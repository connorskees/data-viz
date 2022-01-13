use std::ops::{Add, Div, Mul, Rem, Sub};

use minifb::{Key, Window, WindowOptions};
use num::{One, ToPrimitive, Zero};

pub struct Bitmap {
    pub buffer: Vec<u32>,
    width_px: usize,
    height_px: usize,
}

const RED: u32 = 0xff_ff_00_00;
const WHITE: u32 = 0xff_ff_ff_ff;
const BLACK: u32 = 0x00_00_00_00;

impl Bitmap {
    pub fn new(width_px: usize, height_px: usize) -> Self {
        Self {
            buffer: vec![WHITE; width_px * height_px],
            width_px,
            height_px,
        }
    }

    pub fn draw_vertical_line(&mut self, x: usize) {
        for y in 0..self.height_px {
            let idx = x + y * self.width_px;

            if let Some(pixel) = self.buffer.get_mut(idx) {
                *pixel = BLACK;
            }
        }
    }

    pub fn draw_horizontal_line(&mut self, y: usize) {
        let start = self.width_px * y;

        if let Some(range) = self.buffer.get_mut(start..(start + self.width_px)) {
            range.copy_from_slice(&vec![BLACK; self.width_px]);
        }
    }

    fn set_pixel(&mut self, x: i32, y: i32, color: u32) {
        if x < 0 || y < 0 {
            return;
        }

        let y = self.height_px - y as usize;
        // assert!(x >= 0, "{}", x);
        // assert!(y >= 0, "{}", y);

        let idx = x as usize + y as usize * self.width_px;
        if let Some(pixel) = self.buffer.get_mut(idx) {
            *pixel = color;
        }
    }

    pub fn draw_point(&mut self, origin: (usize, usize), radius: usize) {
        let origin_x = origin.0 as i32;
        let origin_y = origin.1 as i32;

        let mut radius = radius as i32;

        let mut x = -radius;
        let mut y = 0;
        let mut err = 2 - 2 * radius;

        while x < 0 {
            dbg!(x, y, origin_x, origin_y);
            self.set_pixel(origin_x - x, origin_y + y, RED);
            self.set_pixel(origin_x - y, origin_y - x, RED);
            self.set_pixel(origin_x + x, origin_y - y, RED);
            self.set_pixel(origin_x + y, origin_y + x, RED);

            radius = err;

            if radius <= y {
                y += 1;
                err += y * 2 + 1;
            }

            if radius > x || err > y {
                x += 1;
                err += x * 2 + 1;
            }
        }
    }
}

pub struct CartesianGrid<T: Number> {
    x_axis: Axis<T>,
    y_axis: Axis<T>,
    dist_between_cols: f64,
    dist_between_rows: f64,
    pub bitmap: Bitmap,
}

pub trait Number:
    Sized
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Rem<Output = Self>
    + Div<Output = Self>
    + PartialEq
    + One
    + Zero
    + ToPrimitive
    + Copy
{
}

impl<T> Number for T where
    T: Sized
        + Add<Output = Self>
        + Sub<Output = Self>
        + Mul<Output = Self>
        + Rem<Output = Self>
        + Div<Output = Self>
        + PartialEq
        + One
        + Zero
        + ToPrimitive
        + Copy
{
}

pub struct Axis<T: Number> {
    min: T,
    max: T,
    step: T,
}

impl<T: Number> Axis<T> {
    pub fn new(min: T, max: T, step: T) -> Self {
        Self { min, max, step }
    }

    pub fn num_of_lines(&self) -> T {
        (self.max - self.min) / self.step
            + if (self.max - self.min) % self.step != T::zero() {
                T::one()
            } else {
                T::zero()
            }
    }
}

impl<T: Number> CartesianGrid<T> {
    pub fn new(x_axis: Axis<T>, y_axis: Axis<T>, width: usize, height: usize) -> Self {
        let bitmap = Bitmap::new(width, height);

        let num_of_cols = x_axis.num_of_lines().to_usize().unwrap();
        let dist_between_cols = bitmap.width_px as f64 / num_of_cols as f64;

        let num_of_rows = y_axis.num_of_lines().to_usize().unwrap();
        let dist_between_rows = bitmap.height_px as f64 / num_of_rows as f64;

        let mut grid = Self {
            x_axis,
            y_axis,
            dist_between_cols,
            dist_between_rows,
            bitmap,
        };

        grid.write_to_bitmap();

        grid
    }

    pub fn graph_series(&mut self, series: &Series<T>) {
        for &(coord_x, coord_y) in &series.coords {
            let normalized_x = (coord_x - self.x_axis.min).to_f64().unwrap()
                / (self.x_axis.num_of_lines() * self.x_axis.step)
                    .to_f64()
                    .unwrap();
            let normalized_y = (coord_y - self.y_axis.min).to_f64().unwrap()
                / (self.y_axis.num_of_lines() * self.y_axis.step)
                    .to_f64()
                    .unwrap();

            let origin_x = normalized_x * self.bitmap.width_px as f64;
            let origin_y = normalized_y * self.bitmap.height_px as f64;

            // dbg!(idx, origin, coord_x, self.max, self.num_of_cols * self.step);

            // dbg!(idx, coord_x);
            // let pos_x = (coord_x as f64 / (max_x - min_x) as f64) * self.bitmap.width_px as f64;
            // let pos_y = (coord_y as f64 / (max_y - min_y) as f64) * self.bitmap.height_px as f64;

            // dbg!((coord_x as f64 / (max_x - min_x) as f64), coord_x, max_x);

            // self.bitmap.draw_point((pos_x as usize, pos_y as usize), 5);

            // dbg!(pos, coord_x);

            // let x_offset = match self.x_elements.binary_search(coord_x) {
            //     Ok(idx) => idx * self.dist_between_cols,
            //     Err(e) => {}
            // };

            // let y_offset = match self.y_elements.binary_search(coord_y) {
            //     Ok(idx) => idx * self.dist_between_rows,
            //     Err(e) => todo!(),
            // };

            // dbg!(x_offset, y_offset);

            // self.bitmap.draw_point((x_offset, y_offset), 5);
            self.bitmap
                .draw_point((origin_x as usize, origin_y as usize), 3);
        }

        // todo!()
    }

    fn write_to_bitmap(&mut self) {
        let mut pos = self.bitmap.width_px as f64;

        loop {
            dbg!(pos);
            self.bitmap.draw_vertical_line(pos as usize);

            if pos < self.dist_between_cols {
                break;
            }

            pos -= self.dist_between_cols;
        }

        let mut pos = self.bitmap.height_px as f64;

        loop {
            dbg!(pos);
            self.bitmap.draw_horizontal_line(pos as usize);

            if pos < self.dist_between_rows {
                break;
            }

            pos -= self.dist_between_rows;
        }
    }
}

pub struct Series<T: Number> {
    pub coords: Vec<(T, T)>,
}

fn main() {
    let WIDTH = 500;
    let HEIGHT = 500;

    let grid = GraphBuilder::new(500, 500).x_min(5.0).finish();

    // let mut grid = CartesianGrid::new(Axis::new(0, 23, 6), Axis::new(0, 23, 6), WIDTH, HEIGHT);

    // let series = Series {
    //     coords: vec![
    //         (0, 0),
    //         (1, 1),
    //         (2, 2),
    //         (3, 3),
    //         (4, 4),
    //         (5, 5),
    //         (6, 6),
    //         (15, 15),
    //         (6, 2),
    //         (8, 2),
    //         (2, 4),
    //     ],
    // };

    // grid.graph_series(&series);

    let mut window =
        Window::new("Coordinate Plane", WIDTH, HEIGHT, WindowOptions::default()).unwrap();

    window.limit_update_rate(Some(std::time::Duration::from_micros(16600)));

    while window.is_open() && !window.is_key_down(Key::Escape) {
        window
            .update_with_buffer(&grid.bitmap.buffer, WIDTH, HEIGHT)
            .unwrap();
    }
}

struct GraphBuilder<T: Number> {
    x_axis: Axis<T>,
    y_axis: Axis<T>,
    series: Series<T>,
    width: usize,
    height: usize,
}

impl<T: Number> GraphBuilder<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            series: Series { coords: Vec::new() },
            x_axis: Axis::new(T::zero(), T::zero(), T::one()),
            y_axis: Axis::new(T::zero(), T::zero(), T::one()),
            width,
            height,
        }
    }

    // pub fn with_series(series: Series<T>) -> Self<T> {
    //     let x_max = series.max();

    //     Self {
    //         series: Series { coords: Vec::new(), }
    //         x_axis: Axis::new(T::zero(), x_max, T::one()),
    //         y_axis: Axis::new(T::zero(), T::zero(), T::one()),
    //     }
    // }

    fn x_min(mut self, min: T) -> Self {
        self.x_axis.min = min;

        self
    }

    fn finish(self) -> CartesianGrid<T> {
        CartesianGrid::new(self.x_axis, self.y_axis, self.width, self.height)
    }
}
