#include "mandelbrotViewer.h"
#include <string.h>
#include <iostream>
#include <ctime>

//initialize a couple of global objects
sf::Mutex mutex1;
sf::Mutex mutex2;

//Constructor
MandelbrotViewer::MandelbrotViewer(int res) {
    resolution = res;
    height = width = res;

    //create the window and view, then give them to the pointers
    static sf::RenderWindow win(sf::VideoMode(resolution, resolution), "Mandelbrot Explorer");
    static sf::View vw(sf::FloatRect(0, 0, resolution, resolution));
    window = &win;
    view = &vw;

    //initialize the viewport
    view->setViewport(sf::FloatRect(0, 0, 1, 1));
    window->setView(*view);
    
    //cap the framerate
    framerateLimit = 60;
    window->setFramerateLimit(framerateLimit);

    //initialize the mandelbrot parameters
    resetMandelbrot();

    //initialize the image
    texture.create(resolution, resolution);
    image.create(resolution, resolution, sf::Color::Black);
    sprite.setTexture(texture);
    scheme = 1;
    initPalette(); 

    size_t size = resolution;
    std::vector< std::vector<int> > array(size, std::vector<int>(size));
    image_array = array;
}

MandelbrotViewer::~MandelbrotViewer() { }

//Accessors
sf::Vector2i MandelbrotViewer::getMousePosition() {
    return sf::Mouse::getPosition(*window);
}

//return the center of the area of the complex plane
sf::Vector2f MandelbrotViewer::getMandelbrotCenter() {
    sf::Vector2f center;
    center.x = area.left + area.width/2.0;
    center.y = area.top + area.height/2.0;
    return center;
}

//gets the next event from the viewer
bool MandelbrotViewer::getEvent(sf::Event& event) {
    return window->waitEvent(event);
}

//checks if the window is open
bool MandelbrotViewer::isOpen() {
    return window->isOpen();
}

//this is a setter function to change the color scheme
//it also handles all the regeneration and refreshing
void MandelbrotViewer::setColorScheme(int newScheme) {
    scheme = newScheme;
    initPalette();
    changeColor();
    updateMandelbrot();
    refreshWindow();
}

//Functions to change parameters of mandelbrot

//regenerates the image with the new color multiplier, without regenerating
//the mandelbrot. Does not update the image (use updateImage())
void MandelbrotViewer::changeColor() {
    for (int i=0; i<resolution; i++) {
        for (int j=0; j<resolution; j++) {
            image.setPixel(j, i, findColor(image_array[i][j]));
        }
    }
}

//changes the parameters of the mandelbrot: sets new center and zooms accordingly
//does not regenerate or update the image
void MandelbrotViewer::changePos(sf::Vector2<double> new_center, double zoom_factor) {
    area.width = area.width * zoom_factor;
    area.height = area.height * zoom_factor;
    area.left = new_center.x - area.width / 2.0;
    area.top = new_center.y - area.height / 2.0;
    //NOTE: this is a relative zoom
}

//similar to changePos, but it's an absolute zoom and it only changes the view
//instead of setting new parameters to regenerate the mandelbrot
void MandelbrotViewer::changePosView(sf::Vector2f new_center, double zoom_factor) {
    //reset the view so that we can apply an absolute zoom, instead of relative
    resetView();

    //set new center and zoom
    view->setCenter(new_center);
    view->zoom(zoom_factor);

    //reload the new view
    window->setView(*view);
}

// Resize window event
void MandelbrotViewer::resizeWindow(int newHeight, int newWidth) {
    // Getting the center
    sf::Vector2f center = getMandelbrotCenter();

    // Getting the delta width
    int deltaWidth  = newWidth  - width,
        deltaHeight = newHeight - height;
    width = newWidth;
    height = newHeight;

    // Changing the corner data
    sf::Vector2f cornerDelta;
    cornerDelta.x = deltaWidth;
    cornerDelta.y = deltaHeight;
    sf::Vector2<double> complexDelta = pixelToComplex(cornerDelta);
    area.top -= complexDelta.x;
    area.left -= complexDelta.y;
    // Changing the size data
    sf::Vector2f sizeDelta;
    sizeDelta.x = width;
    sizeDelta.y = height;
    complexDelta = pixelToComplex(sizeDelta);
    area.width = complexDelta.x;
    area.height= complexDelta.y;
}

//generate the mandelbrot
void MandelbrotViewer::generate() {

    //make sure it starts at line 0
    nextLine = 0;

    sf::Thread thread1(&MandelbrotViewer::genLine, this);
    sf::Thread thread2(&MandelbrotViewer::genLine, this);
    sf::Thread thread3(&MandelbrotViewer::genLine, this);
    sf::Thread thread4(&MandelbrotViewer::genLine, this);

    thread1.launch();
    thread2.launch();
    thread3.launch();
    thread4.launch();

    thread1.wait();
    thread2.wait();
    thread3.wait();
    thread4.wait();

    last_max_iter = max_iter;
}

//this is a private worker thread function. Each thread picks the next ungenerated
//row of pixels, generates it, then starts the next one
void MandelbrotViewer::genLine() {

    int iter, row, column;
    double x, y;
    double x_inc = interpolate(area.width, resolution);
    double y_inc = interpolate(area.height, resolution);
    sf::Color color;

    while(true) {

        //the mutex avoids multiple threads writing to variables at the same time,
        //which can corrupt the data
        mutex1.lock();
        row = nextLine++; //get the next ungenerated line
        mutex1.unlock();

        //if all the rows have been generated, stop it from generating outside the bounds
        //of the image
        if (row >= resolution) return;

        //calculate the row height in the complex plane
        y = area.top + row * y_inc;

        //now loop through and generate all the pixels in that row
        for (column = 0; column < resolution; column++) {

            // Check if we increased iterations and if the pixel already diverged
            if ( last_max_iter < max_iter && image_array[row][column] < last_max_iter ) {
                iter = image_array[row][column];
            } // Check if we decreased iterations and if the pixel already converged
            else if ( last_max_iter > max_iter && image_array[row][column] > max_iter) {
                iter = image_array[row][column];
            } // Check if we zoomed, or didn't change iterations which means we need to recalculate the whole thing
            else {
                //calculate the next x coordinate of the complex plane
                x = area.left + column * x_inc;
                iter = escape(x, y);
            }

            //mutex this too so that the image is not accessed multiple times simultaneously
            mutex2.lock();
            image.setPixel(column, row, findColor(iter));
            image_array[row][column] = iter;
            mutex2.unlock();
        }
    }
}

//Reset/update functions:

//resets the mandelbrot to generate the starting area
void MandelbrotViewer::resetMandelbrot() {
    area.left = -1.5;
    area.top = -1.0;
    area.width = 2;
    area.height = 2;
    max_iter = 100;
    last_max_iter = 100;
    color_multiple = 1;
}

//refreshes the window: clear, draw, display
void MandelbrotViewer::refreshWindow() {
    window->clear(sf::Color::Black);
    window->setView(*view);
    window->draw(sprite);
    window->display();
}

//reset the view to display the entire image
void MandelbrotViewer::resetView() {
    view->reset(sf::FloatRect(0, 0, resolution, resolution));
    window->setView(*view);
}

//close the window
void MandelbrotViewer::close() {
    window->close();
}

//update the mandelbrot image (use the already generated image to update the
//texture, so the next time the screen updates it will be displayed
void MandelbrotViewer::updateMandelbrot() {
    texture.update(image);
}

void MandelbrotViewer::setWindowActive(bool setting) {
    window->setActive(setting);
}

//saves the currently displayed image to a png with a timestamp in the title
void MandelbrotViewer::saveImage() {
    //set up the timestamp filename
    time_t currentTime = time(0);
    tm* currentDate = localtime(&currentTime);
    char filename[80];
    strftime(filename,80,"%Y-%m-%d.%H-%M-%S",currentDate);
    strcat(filename, ".png");

    //save the image and print confirmation
    image.saveToFile(filename);
    std::cout << "Saved image to " << filename << std::endl;
}

//Converts a vector from pixel coordinates to the corresponding
//coordinates on the complex plane
sf::Vector2<double> MandelbrotViewer::pixelToComplex(sf::Vector2f pix) {
    sf::Vector2<double> comp;
    comp.x = area.left + pix.x * interpolate(area.width, resolution);
    comp.y = area.top + pix.y * interpolate(area.height, resolution);
    return comp;
}

//this function calculates the escape-time of the given coordinate
//it is the brain of the mandelbrot program: it does the work to
//make the pretty pictures :)
int MandelbrotViewer::escape(double x0, double y0) {
    double x = 0, y = 0, x_check = 0, y_check = 0;
    int iter = 0, period = 2;

    double x_square = 0;
    double y_square = 0;

    //this is a specialized version of z = z^2 + c. It only does three multiplications,
    //instead of the normal six. Multplications are very costly with such high precision
    while(period < max_iter) {
        x_check = x;
        y_check = y;
        period += period;

        if (period > max_iter) period = max_iter;
        for (; iter < period; iter++) {
            y = x * y;
            y += y; //multiply by two
            y += y0;
            x = x_square - y_square + x0;

            x_square = x*x;
            y_square = y*y;

            //if the magnitued is greater than 2, it will escape
            if (x_square + y_square > 4.0) return iter;

            //another optimization: it checks if the new 'z' is a repeat. If so,
            //it knows that it is in a loop and will not escape
            if ((x == x_check) && (y == y_check)){
                return max_iter;
            }
        }
    }
    return max_iter;
}

//findColor uses the number of iterations passed to it to look up a color in the palette
sf::Color MandelbrotViewer::findColor(int iter) {
    int i = fmod(iter * color_multiple, 255);
    sf::Color color;
    if (iter >= max_iter) color = sf::Color::Black;
    else {
        color.r = palette[0][i];
        color.g = palette[1][i];
        color.b = palette[2][i];
    }
    return color;
}

//This is for initPalette, it makes sure the given number is between 0 and 255
int coerce(int number) {
    if (number > 255) number = 255;
    else if (number < 0) number = 0;
    return number;
}

//Sets up the palette array
void MandelbrotViewer::initPalette() {
    //define some non-standard colors
    sf::Color orange;
    orange.r = 255;
    orange.g = 165;
    orange.b = 0;
    sf::Color pink;
    pink.r = 255;
    pink.g = 135;
    pink.b = 135;
    sf::Color dark_red;
    dark_red.r = 209;
    dark_red.g = 2;
    dark_red.b = 2;
    sf::Color dark_green;
    dark_green.r = 8;
    dark_green.g = 172;
    dark_green.b = 8;
    sf::Color light_green;
    light_green.r = 115;
    light_green.g = 255;
    light_green.b = 115;
    switch (scheme) {
        //scheme one is black:blue:white:orange:black
        case 1:
            smoosh(sf::Color::Black, sf::Color::Blue, 0, 64);
            smoosh(sf::Color::Blue, sf::Color::White, 64, 144);
            smoosh(sf::Color::White, orange, 144, 196);
            smoosh(orange, sf::Color::Black, 196, 256);
            break;
        //scheme two is black:green:blue:black
        case 2:
            smoosh(sf::Color::Black, sf::Color::Green, 0, 85);
            smoosh(sf::Color::Green, sf::Color::Blue, 85, 170);
            smoosh(sf::Color::Blue, sf::Color::Black, 170, 256);
            break;
        //scheme three is black:red:orange:black
        case 3:
            smoosh(sf::Color::Black, sf::Color::Red, 0, 180);
            smoosh(sf::Color::Red, orange, 180, 215);
            smoosh(orange, sf::Color::Black, 215, 256);
            break;
        //scheme four is black:cyan:white:black
        case 4:
            smoosh(sf::Color::Black, sf::Color::Cyan, 0, 110);
            smoosh(sf::Color::Cyan, sf::Color::White, 110, 220);
            smoosh(sf::Color::White, sf::Color::Black, 220, 256);
            break;
        //scheme five is red:orange:yellow:green:blue:magenta:red
        case 5:
            smoosh(sf::Color::Red, orange, 0, 42);
            smoosh(orange, sf::Color::Yellow, 42, 84);
            smoosh(sf::Color::Yellow, sf::Color::Green, 84, 127);
            smoosh(sf::Color::Green, sf::Color::Blue, 127, 170);
            smoosh(sf::Color::Blue, sf::Color::Magenta, 170, 212);
            smoosh(sf::Color::Magenta, sf::Color::Red, 212, 256);
            break;
        //scheme six is dark_red:pink, then dark_green:light_green
        case 6:
            smoosh(dark_red, pink, 0, 128);
            smoosh(dark_green, light_green, 128, 256);
            break;
        //scheme seven is smoky... black:white
        case 7:
            smoosh(sf::Color::White, sf::Color::Black, 0, 256);
            break;
        default:
            int r, g, b;
            for (int i = 0; i <= 255; i++) {
                r = (int) (23.45 - 1.880*i + 0.0461*i*i - 0.000152*i*i*i);
                g = (int) (17.30 - 0.417*i + 0.0273*i*i - 0.000101*i*i*i);
                b = (int) (25.22 + 7.902*i - 0.0681*i*i + 0.000145*i*i*i);

                palette[0][i] = coerce(r);
                palette[1][i] = coerce(g);
                palette[2][i] = coerce(b);
            }
    }
}

//Smooshes two colors together, and writes them to the palette in the specified range
void MandelbrotViewer::smoosh(sf::Color c1, sf::Color c2, int min, int max) {
    int range = max-min;
    float r_inc = interpolate(c1.r, c2.r, range);
    float g_inc = interpolate(c1.g, c2.g, range);
    float b_inc = interpolate(c1.b, c2.b, range);

    //loop through the palette setting new colors
    for (int i=0; i < range; i++) {
        palette[0][min+i] = (int) (c1.r + i * r_inc);
        palette[1][min+i] = (int) (c1.g + i * g_inc);
        palette[2][min+i] = (int) (c1.b + i * b_inc);
    }
}
