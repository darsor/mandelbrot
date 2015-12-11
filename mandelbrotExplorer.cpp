#include "mandelbrotViewer.h"

//global stuff
int iterations = 100;
int resolution;

//this struct holds the parameters for the zoom function,
//since it needs to be threaded and thus cannot accept arguments
struct zoomParameters {
    MandelbrotViewer *viewer;
    sf::Vector2f oldc;
    sf::Vector2f newc;
    double zoom;
    int frames;
} param;

//returns range increments to get from min to max
double interpolate(double min, double max, int range) { return (max-min)/range; }

//function prototypes
void handleKeyboard(MandelbrotViewer *brot, sf::Event *event);
void handleZoom(MandelbrotViewer *brot, sf::Event *event);
void handleDrag(MandelbrotViewer *brot, sf::Event *event);
void zoom();

int main() {
    resolution = 820;
    
    //create the mandelbrotviewer instance
    MandelbrotViewer brot(resolution);

    //initialize the image
    brot.resetMandelbrot();
    brot.generate();
    brot.updateMandelbrot();
    brot.refreshWindow();

    //point the zoom function to the 'brot' instance
    param.viewer = &brot;

    //create an event to test for input
    sf::Event event;

    //main window loop
    while (brot.isOpen()) {

        brot.getEvent(event);

        //this big switch statement handles all types of input
        switch (event.type) {

            //if the event is a keypress
            case sf::Event::KeyPressed:
                handleKeyboard(&brot, &event);
                break;

            //if the event is a mousewheel scroll, zoom
            case sf::Event::MouseWheelScrolled:
                handleZoom(&brot, &event);
                break;

            //if the event is a click, drag the view:
            case sf::Event::MouseButtonPressed:
                handleDrag(&brot, &event);
                break;

            default:
                break;
        } //end event switch

    } //end main window loop

    return 0;
} //end main

//this function handles all keyboard input
void handleKeyboard(MandelbrotViewer *brot, sf::Event *event) {
    float color_inc = interpolate(0, 1, 30);
    switch(event->key.code) {
        //if Q, quit and close the window
        case sf::Keyboard::Q:
            brot->close();
            break;
        //if up arrow, increase iterations
        case sf::Keyboard::Up:
            iterations += 30;
            brot->setIterations(iterations);
            brot->generate();
            brot->updateMandelbrot();
            brot->refreshWindow();
            break;
        //if down arrow, decrease iterations
        case sf::Keyboard::Down:
            iterations -= 30;
            if (iterations < 100) iterations = 100;
            brot->setIterations(iterations);
            brot->generate();
            brot->updateMandelbrot();
            brot->refreshWindow();
            break;
        //if right arrow, increase color_multiple until released
        case sf::Keyboard::Right:
            color_inc = interpolate(0, 1, 25);
            while (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) {
                brot->setColorMultiple(brot->getColorMultiple() + color_inc);
                brot->changeColor();
                brot->updateMandelbrot();
                brot->refreshWindow();
            }
            break;
        //if left arrow, decrease color_multiple until released
        case sf::Keyboard::Left:
            color_inc = interpolate(1, 0, 25);
            while (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) {
                if (brot->getColorMultiple() > 1) {
                    brot->setColorMultiple(brot->getColorMultiple() + color_inc);
                    brot->changeColor();
                    brot->updateMandelbrot();
                    brot->refreshWindow();
                }
            }
            break;
        //if it's a 1, change to color scheme 1
        case sf::Keyboard::Num1:
            brot->setColorScheme(1);
            break;
        //if it's a 2, change to color scheme 2
        case sf::Keyboard::Num2:
            brot->setColorScheme(2);
            break;
        //if it's a 3, change to color scheme 3
        case sf::Keyboard::Num3:
            brot->setColorScheme(3);
            break;
        //if it's a 4, change to color scheme 4
        case sf::Keyboard::Num4:
            brot->setColorScheme(4);
            break;
        //if it's a 5, change to color scheme 5
        case sf::Keyboard::Num5:
            brot->setColorScheme(5);
            break;
        //if it's a 6, change to color scheme 6
        case sf::Keyboard::Num6:
            brot->setColorScheme(6);
            break;
        //if it's a 7, change to color scheme 7
        case sf::Keyboard::Num7:
            brot->setColorScheme(7);
            break;
        //if R, reset the mandelbrot to the starting image
        case sf::Keyboard::R:
            brot->resetMandelbrot();
            brot->resetView();
            brot->generate();
            brot->updateMandelbrot();
            brot->refreshWindow();
            break;
        //if S, save the current image
        case sf::Keyboard::S:
            brot->saveImage();
            break;
        //end the keypress switch
        default:
            break;
    }
}

//this function handles scroll wheel input
void handleZoom(MandelbrotViewer *brot, sf::Event *event){

    //get some vectors ready to calculate the zoom
    sf::Vector2f old_center;
    sf::Vector2f new_center;

    //zoom in with the mouse coordinates as the new center
    //set the zoom frames up so that it is more animated
    param.frames = 30;
    old_center = brot->getViewCenter();
    new_center.x = event->mouseWheelScroll.x;
    new_center.y = event->mouseWheelScroll.y;

    //if it's an upward scroll, get ready to zoom in
    if (event->mouseWheelScroll.delta > 0) {
        brot->changePos(brot->pixelToComplex(new_center), 0.5);
        param.zoom = 0.5;
    } //if it's a downward scroll, get ready to zoom out
    else if (event->mouseWheelScroll.delta < 0) {
        brot->changePos(brot->pixelToComplex(new_center), 2.0);
        param.zoom = 2.0;
    }

    //set the parameters for the zoom
    param.oldc = old_center;
    param.newc = new_center;

    //SFML doesn't like having a window active in more than one thread
    brot->setWindowActive(false);

    //start zooming with a worker thread, so that it can generate
    //the new image while it's zooming
    sf::Thread thread(&zoom);
    thread.launch();

    //start generating while it's zooming
    brot->generate();

    //wait for the thread to finish (wait for the zoom to finish)
    thread.wait();

    //now reactivate the window
    brot->setWindowActive(true);

    //now display the new mandelbrot
    brot->updateMandelbrot();
    brot->resetView();
    brot->refreshWindow();
}

void handleDrag(MandelbrotViewer *brot, sf::Event *event) {
    
    int framerateLimit = brot->getFramerate();

    //get some vectors ready for calculations
    sf::Vector2f old_center;
    sf::Vector2f new_center;
    sf::Vector2f old_position;
    sf::Vector2f new_position;
    sf::Vector2f difference;
    sf::Vector2i temp;

    //set drag frames low so that it will move in real time
    param.frames = 2;

    //save the old center/mouse_position as reference
    temp = brot->getMousePosition();
    old_position.x = temp.x;
    old_position.y = temp.y;
    old_center = brot->getViewCenter();
    param.oldc.x = old_center.x;
    param.oldc.y = old_center.y;

    //set zoom to 1 so that it doesn't change, only drags
    param.zoom = 1.0;

    //set the framrate very high so that it will drag in real time
    brot->setFramerate(500);

    while (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {

        //get the new mouse position
        temp = brot->getMousePosition();
        new_position.x = temp.x;
        new_position.y = temp.y;

        //get the difference between the old and new positions
        difference.x = new_position.x - old_position.x;
        difference.y = new_position.y - old_position.y;

        //calculate the new center
        new_center = old_center - difference;
        param.newc.x = new_center.x;
        param.newc.y = new_center.y;

        //run the zoom function, but since param.zoom is set to 1,
        //it will only drag
        zoom();

        //start over again until the mouse is released
        param.oldc = param.newc;
    }

    //set the framerate back to default
    brot->setFramerate(framerateLimit);

    //now regenerate the mandelbrot at the new position.
    //This can't be threaded like zoom, because it doesn't
    //know where to generate at until it's done dragging
    temp = brot->getMousePosition();
    new_position.x = temp.x;
    new_position.y = temp.y;

    //get the difference of how far it's moved
    difference.x = new_position.x - old_position.x;
    difference.y = new_position.y - old_position.y;

    //set the old center
    old_center = sf::Vector2f(resolution/2, resolution/2);

    //calculate the new center
    new_center = old_center - difference;

    brot->changePos(brot->pixelToComplex(new_center), 1.0);
    brot->generate();
    brot->resetView();
    brot->updateMandelbrot();
    brot->refreshWindow();
}

//animates the zoom of the viewer (so the viewer zooms while the mandelbrot is generating)
//uses the struct param as parameters
void zoom() {
    param.viewer->setWindowActive(true);
    double inc_drag_x = interpolate(param.oldc.x, param.newc.x, param.frames);
    double inc_drag_y = interpolate(param.oldc.y, param.newc.y, param.frames);
    double inc_zoom = interpolate(1.0, param.zoom, param.frames);

    //animate the zoom
    for (int i=0; i<param.frames; i++) {
        param.newc.x = param.oldc.x + i * inc_drag_x;
        param.newc.y = param.oldc.y + i * inc_drag_y;
        param.viewer->changePosView(param.newc, 1 + i * inc_zoom);
        param.viewer->refreshWindow();
    }
    param.viewer->setWindowActive(false);
}
