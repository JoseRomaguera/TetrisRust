extern crate piston_window;
use piston_window::*;
use std::time::{SystemTime, UNIX_EPOCH};

const TICK_DELAY: f32 = 0.5;
const TICK_RATE: f32 = 1.0 / TICK_DELAY;
const ROWS: usize = 23;
const COLUMNS: usize = 10;
const DEFEAT_HEIGHT: i32 = 20;

const CELL_SIZE: f32 = 30.0;
const GRID_WIDTH: f32 = CELL_SIZE * COLUMNS as f32;
const GRID_HEIGHT: f32 = CELL_SIZE * ROWS as f32;

const CELL_VISUAL_SCALE: f32 = 0.9;
const BORDER_THICKNESS: f32 = 10.0;
const BORDER_COLOR: [f32; 4] = [0.0, 0.8, 1.0, 1.0];
const DEFEAT_COLOR: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum CellValue { None, Red, Green, Blue, count }

const CELL_COLORS: &[[f32; 4]] = &[[0.1, 0.1, 0.1, 1.0], [1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0]];

const CELL_SHAPES: &[[u8; 3 * 3]] = &[
    [
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
    ],
    [
        0, 0, 0,
        1, 1, 0,
        0, 1, 1,
    ],
    [
        0, 0, 0,
        0, 1, 1,
        1, 1, 0,
    ],
    [
        0, 0, 0,
        0, 0, 1,
        1, 1, 1,
    ],
    [
        0, 0, 0,
        1, 0, 0,
        1, 1, 1,
    ],
    [
        0, 0, 0,
        1, 1, 0,
        1, 1, 0,
    ],
];

struct Game {
    screen_width: i32,
    screen_height: i32,
    dt: f32,
    tick_timer: f32,
    tick: i32,
    speedup: bool,

    score: i32,
    finished: bool,

    seed: u64,
    seed_mask: u64,

    grid: [CellValue; ROWS * COLUMNS],
    shape_index: i32,
    shape_value: CellValue,
    shape_position: (i32, i32),
    shape_rotation: i32,
}

impl Game {
    fn new() -> Self {
        Game {
            screen_width: 1,
            screen_height: 1,
            dt: 0.001,
            tick_timer: 0.0,
            tick: 0,
            speedup: false,
            score: 0,
            finished: false,
            seed: 1,
            seed_mask: 1,
            grid: [CellValue::None; ROWS * COLUMNS],
            shape_index: 0,
            shape_value: CellValue::Red,
            shape_position: (0, 0),
            shape_rotation: 0,
        }
    }
}

fn get_current_time_seed() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos() as u64
}

fn get_cell_color(value: CellValue) -> [f32; 4] {
    let index = value as i32;
    if index < 0 || index >= CELL_COLORS.len() as i32 { return [0.5, 0.5, 0.5, 1.0]; }
    return CELL_COLORS[index as usize];
}

fn cell_value_from_index(value: u32) -> CellValue {
    match value {
        0 => CellValue::None,
        1 => CellValue::Red,
        2 => CellValue::Blue,
        3 => CellValue::Green,
        _ => { assert!(false, "Unknown cell value"); CellValue::None }
    }
}

fn calculate_cell_center(col: i32, row: i32) -> (f32, f32) {
    let x = col as f32 * CELL_SIZE + CELL_SIZE * 0.5;
    let y = row as f32 * CELL_SIZE + CELL_SIZE * 0.5;
    return (x, y);
}

// From: https://github.com/svaarala/duktape/blob/master/misc/splitmix64.c
fn random_splitmix64(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *seed;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

fn u64_random(game: &mut Game) -> u64 {
    let mut seed: u64 = game.seed ^ game.seed_mask;
    let result: u64 = random_splitmix64(&mut seed);
    game.seed_mask = seed;
    return result;
}

fn u32_random(game: &mut Game) -> u32 {
	return u64_random(game) as u32;
}
fn u32_random_max(game: &mut Game, max: u32) -> u32 {
	return u32_random(game) % max;
}
fn u32_random_range(game: &mut Game, min: u32, max: u32) -> u32 {
	return min + (u32_random(game) % (((max as i32) - (min as i32)) as u32));
}
fn f32_random(game: &mut Game) -> f32 {
	let v: u32 = u32_random(game) & 0xFFFFFF;
	return (v as f32) / (0xFFFFFF as f32);
}
fn f32_random_max(game: &mut Game, max: f32) -> f32 {
	return f32_random(game) * max;
}
fn f32_random_range(game: &mut Game, min: f32, max: f32) -> f32 {
	return min + f32_random(game) * (max - min);
}

fn draw_quad(game: &Game, x: f32, y: f32, width: f32, height: f32, color: [f32; 4], c: Context, g: &mut G2d)
{
    let mut x0 = x;
    let mut y0 = y;

    // Center grid
    x0 -= GRID_WIDTH * 0.5;
    y0 -= GRID_HEIGHT * 0.5;

    // Change coord system: From mine to Pistons
    let mut x0 = x0 + game.screen_width as f32 * 0.5 - width * 0.5;
    let mut y0 = -y0 + game.screen_height as f32 * 0.5 - height * 0.5;

    rectangle(color, [x0 as f64, y0 as f64, width as f64, height as f64], c.transform, g);
}

fn draw_cell(game: &Game, x: f32, y: f32, color: [f32; 4], c: Context, g: &mut G2d)
{
    draw_quad(game, x, y, CELL_SIZE * CELL_VISUAL_SCALE, CELL_SIZE * CELL_VISUAL_SCALE, color, c, g);
}

fn draw_grid(game: &Game, c: Context, g: &mut G2d)
{
    // Border
    {
        let th = BORDER_THICKNESS;

        draw_quad(game, -th * 0.5, GRID_HEIGHT * 0.5, th, GRID_HEIGHT, BORDER_COLOR, c, g);
        draw_quad(game, GRID_WIDTH + th * 0.5, GRID_HEIGHT * 0.5, th, GRID_HEIGHT, BORDER_COLOR, c, g);

        draw_quad(game, GRID_WIDTH * 0.5, -th * 0.5, GRID_WIDTH + th * 2.0, th, BORDER_COLOR, c, g);
        draw_quad(game, GRID_WIDTH * 0.5, GRID_HEIGHT + th * 0.5, GRID_WIDTH + th * 2.0, th, BORDER_COLOR, c, g);
    }

    for row in 0..ROWS {
        for col in 0..COLUMNS {
            let v = game.grid[col + row * COLUMNS];
            let color = get_cell_color(v);

            let pos = calculate_cell_center(col as i32, row as i32);
            draw_cell(game, pos.0, pos.1, color, c, g);
        }
    }

    // Defeat Line
    draw_quad(game, GRID_WIDTH * 0.5, DEFEAT_HEIGHT as f32 * CELL_SIZE, GRID_WIDTH, 5.0, DEFEAT_COLOR, c, g);
}

fn calculate_shape_cell(col0: i32, row0: i32, rotation: i32) -> (i32, i32) {
    let mut col = col0;
    let mut row = row0;

    match rotation {
        0 => {},
        1 => { let aux = col; col = row; row = 2 - aux; },
        2 => { row = 2 - row; col = 2 - col; },
        3 => { let aux = col; col = 2 - row; row = aux; },
        _ => {}
    }
    return (col, row);
}

fn draw_shape(game: &Game, shape: [u8; 3 * 3], col: i32, row: i32, rotation: i32, value: CellValue, c: Context, g: &mut G2d)
{
    let color = get_cell_color(value);

    for row0 in 0..3 {
        for col0 in 0..3 {
            let filled = shape[(col0 + (2 - row0) * 3) as usize] != 0;
            if !filled {continue;}
            let relative_cell = calculate_shape_cell(col0, row0, rotation);
            let pos = calculate_cell_center(col + relative_cell.0 as i32, row + relative_cell.1 as i32);
            draw_cell(game, pos.0, pos.1, color, c, g);
        }
    }
}

fn assign_new_shape(game: &mut Game)
{
    game.shape_index = u32_random_max(game, CELL_SHAPES.len() as u32) as i32;
    game.shape_value = cell_value_from_index(u32_random_range(game, 1, CellValue::count as u32));
    game.shape_position = (COLUMNS as i32 / 2 - 1, ROWS as i32 - 3);
    game.shape_rotation = 0;
}

fn set_cell(game: &mut Game, col: i32, row: i32, value: CellValue)
{
    if col < 0 || col >= COLUMNS as i32 || row < 0 || row >= ROWS as i32 { 
        assert!(false, "Cell out of bounds");
        return;
    }
    
    game.grid[(col + row * COLUMNS as i32) as usize] = value;
}

fn is_row_completed(game: &mut Game, row: i32) -> bool {
    for col in 0..COLUMNS {
        let v = game.grid[col + (row as usize) * COLUMNS];
        if v == CellValue::None { return false }
    }
    true
}

fn erase_completed_rows(game: &mut Game)
{
    let mut row = 0;
    while row < ROWS
    {
        if is_row_completed(game, row as i32) {
            for row0 in (row + 1)..ROWS {
                for col in 0..COLUMNS {
                    let dst_index = col + (row0 - 1) * COLUMNS;
                    let src_index = col + (row0 + 0) * COLUMNS;
                    game.grid[dst_index] = game.grid[src_index];
                }
            }
            game.score += 1;
        }
        else {
            row += 1;
        }
    }
}

fn apply_shape(game: &mut Game)
{
    let shape = CELL_SHAPES[game.shape_index as usize];

    for row0 in 0..3 {
        for col0 in 0..3 {
            let filled = shape[(col0 + (2 - row0) * 3) as usize] != 0;
            if !filled {continue;}

            let relative_cell = calculate_shape_cell(col0, row0, game.shape_rotation);
            let col = game.shape_position.0 + relative_cell.0 as i32;
            let row = game.shape_position.1 + relative_cell.1 as i32;

            if row >= DEFEAT_HEIGHT { game.finished = true; }

            set_cell(game, col, row, game.shape_value);
        }
    }

    erase_completed_rows(game);
    assign_new_shape(game);
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum CollideResult { None, HorizontalLimit, Bottom, Cell, Up }

fn cell_collides(game: &mut Game, col: i32, row: i32) -> CollideResult
{
    if col < 0 || col >= COLUMNS as i32 { return CollideResult::HorizontalLimit; }
    if row < 0 { return CollideResult::Bottom; }
    if row >= ROWS as i32 { return CollideResult::Up; }
    if game.grid[col as usize + row as usize * COLUMNS] != CellValue::None { return CollideResult::Cell; }
    return CollideResult::None;
}

fn shape_collides(game: &mut Game, shape: [u8; 3 * 3], col: i32, row: i32, rotation: i32) -> CollideResult
{
    let mut collide = CollideResult::None;
    for row0 in 0..3 {
        for col0 in 0..3 {
            let filled = shape[(col0 + (2 - row0) * 3) as usize] != 0;
            if !filled {continue;}
            let relative_cell = calculate_shape_cell(col0, row0, rotation);
            let collide0 = cell_collides(game, col + relative_cell.0 as i32, row + relative_cell.1 as i32);
            if collide0 > collide { collide = collide0; }
        }
    }
    return collide;
}

fn advance_tick(game: &mut Game)
{
    let result = shape_collides(game, CELL_SHAPES[game.shape_index as usize], game.shape_position.0, game.shape_position.1 - 1, game.shape_rotation);
    if result != CollideResult::None && result != CollideResult::HorizontalLimit {
        apply_shape(game);
    }
    else {
        game.shape_position.1 -= 1;
    }
    game.tick += 1;
}

fn rotate_shape(game: &mut Game) {
    game.shape_rotation = (game.shape_rotation + 1) % 4;

    // Handle collisions post rotation
    let dir: i32 = if game.shape_position.0 > (COLUMNS as i32 / 2) {-1} else {1};
    while true {
        let result = shape_collides(game, CELL_SHAPES[game.shape_index as usize], game.shape_position.0, game.shape_position.1, game.shape_rotation);
        match result {
            CollideResult::HorizontalLimit => { game.shape_position.0 += dir; }
            CollideResult::None => { break; }
            _ => { apply_shape(game); break; }
        }
    }
}

fn move_shape(game: &mut Game, dir: i32) {
    let result = shape_collides(game, CELL_SHAPES[game.shape_index as usize], game.shape_position.0 + dir, game.shape_position.1, game.shape_rotation);

    match result {
        CollideResult::HorizontalLimit => { }
        CollideResult::None => { game.shape_position.0 += dir; }
        _ => { apply_shape(game); }
    }
}

fn reset(game: &mut Game) {
    game.tick_timer = 0.0;
    game.tick = 0;
    game.speedup = false;

    game.score = 0;
    game.finished = false;

    game.seed = get_current_time_seed();
    game.seed_mask = 1;

    for i in 0..ROWS * COLUMNS {
        game.grid[i] = CellValue::None;
    }
    assign_new_shape(game);
}

fn main() {
    let mut window: PistonWindow = WindowSettings::new("Tetris", [1000, 1000])
        .exit_on_esc(true)
        .build()
        .unwrap();
    
    let mut game = Game::new();

    reset(&mut game);
    
    assign_new_shape(&mut game);

    while let Some(e) = window.next() {

        if let Some(u) = e.update_args() {
            game.dt = u.dt as f32;   
        }
        else {
            game.dt = 0.001;
        }

        if game.finished { reset(&mut game); }

        if let Some(Button::Keyboard(key)) = e.press_args() {
            if key == Key::W || key == Key::Up { rotate_shape(&mut game); }
            if key == Key::A || key == Key::Left { move_shape(&mut game, -1); }
            if key == Key::D || key == Key::Right { move_shape(&mut game, 1); }
            if key == Key::S || key == Key::Down { game.speedup = true; }
        }

        if let Some(Button::Keyboard(key)) = e.release_args() {
            if key == Key::S || key == Key::Down { game.speedup = false; }
        }

        let mut tick_mult: f32 = 1.0;
        if (game.speedup) { tick_mult *= 10.0; }
        game.tick_timer += game.dt * tick_mult;
        while game.tick_timer > TICK_DELAY {
            game.tick_timer -= TICK_DELAY;
            advance_tick(&mut game);
            game.tick += 1;
        }

        window.draw_2d(&e, |c, g, _| {
            if let Some(vp) = c.viewport {
                game.screen_width = vp.draw_size[0] as i32;
                game.screen_height = vp.draw_size[1] as i32;
            }
            else {
                game.screen_width = 1;
                game.screen_height = 1;
            }

            clear([0.0, 0.0, 0.0, 1.0], g);
            draw_grid(&game, c, g);

            draw_shape(&game, CELL_SHAPES[game.shape_index as usize], game.shape_position.0, game.shape_position.1, game.shape_rotation, game.shape_value, c, g);
        });
    }
}
