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


//generate the mandelbrot
void MandelbrotViewer::generate() {

    //make sure it starts at line 0
    nextLine = 0;

    // create the thread pool
    std::vector<sf::Thread> threadPool;
    for (int i=0; i<THREAD_COUNT; i++) {
        sf::Thread temp_thread(&MandelbrotViewer::genLine, this);
        threadPool.push_back(temp_thread);
    }
    // Launch all the threads
    for (int i=0; i<THREAD_COUNT; i++) {
        threadPool[i].launch();
    }
    // Close the threads
    for (int i=0; i<THREAD_COUNT; i++) {
        threadPool[i].wait();
    }
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

            
            //check if we already know that that point escapes.
            //if it's regenerating after a max_iter change, this saves
            //a lot of time. It's disabled for now (TODO)
            //if (escape_array[row][column] == false) {
            //if (image_array[row][column] != max_iter) {

                //calculate the next x coordinate of the complex plane
                x = area.left + column * x_inc;
                iter = escape(x, y);

                //mutex this too so that the image is not accessed multiple times simultaneously
                mutex2.lock();
                image.setPixel(column, row, findColor(iter));
                image_array[row][column] = iter;
                //if (iter < MAX_ITER) escape_array[row][column] = true;
                mutex2.unlock();
            //} 
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
    //scheme one is black:blue:white:orange:black
    if (scheme == 1) {
        sf::Color orange;
        orange.r = 255;
        orange.g = 165;
        orange.b = 0;
        smoosh(sf::Color::Black, sf::Color::Blue, 0, 64);
        smoosh(sf::Color::Blue, sf::Color::White, 64, 144);
        smoosh(sf::Color::White, orange, 144, 196);
        smoosh(orange, sf::Color::Black, 196, 256);
    } else if (scheme == 2) {
        smoosh(sf::Color::Black, sf::Color::Green, 0, 85);
        smoosh(sf::Color::Green, sf::Color::Blue, 85, 170);
        smoosh(sf::Color::Blue, sf::Color::Black, 170, 256);
    } else if (scheme == 3) {
        smoosh(sf::Color::Black, sf::Color::Red, 0, 200);
        smoosh(sf::Color::Red, sf::Color::Black, 200, 256);
    } else {
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
