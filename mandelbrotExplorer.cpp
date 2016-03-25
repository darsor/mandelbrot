#include "mandelbrotViewer.h"
#include <chrono>
#include <iostream>
#include <thread>

# define PI 3.14159265358979323846

//this struct holds the parameters for the zoom function,
//since it needs to be threaded and thus cannot accept arguments
struct zoomParameters {
    MandelbrotViewer *brot;
    sf::Vector2f oldc;
    sf::Vector2f newc;
    sf::Event event;
    double zoom;
    int frames;
    bool done;
} param;

//returns range increments to get from min to max
double interpolate(double min, double max, int range) { return (max-min)/range; }

//function prototypes
void handleEvent();
void handleKeyboard(MandelbrotViewer *brot, sf::Event *event);
void handleZoom(MandelbrotViewer *brot, sf::Event *event);
void handleDrag(MandelbrotViewer *brot, sf::Event *event);
void handleResize(MandelbrotViewer *brot, sf::Event *event);
double handleRotate();
void handleGenerate();
void eventPoll();
void zoom();

int main() {
    
    //create the mandelbrotviewer instance
    MandelbrotViewer brot(820, 820);

    //initialize the image
    brot.resetMandelbrot();
    brot.generate();
    brot.updateMandelbrot();
    brot.refreshWindow();

    //point the zoom function to the 'brot' instance
    param.brot = &brot;
    param.done = true;

    //create an event to test for input
    sf::Event event;

    //main window loop
    while (brot.isOpen()) {

        brot.waitEvent(param.event);

        handleEvent();

    } //end main window loop

    return 0;
} //end main

//this function handles all events
void handleEvent() {
    //this switch statement handles all types of input
    switch (param.event.type) {

        //if the event is a keypress
        case sf::Event::KeyPressed:
            handleKeyboard(param.brot, &param.event);
            break;

            //if the event is a mousewheel scroll, zoom
        case sf::Event::MouseWheelScrolled:
            handleZoom(param.brot, &param.event);
            break;

            //if the event is a click, drag the view:
        case sf::Event::MouseButtonPressed:
            handleDrag(param.brot, &param.event);
            break;

            //if the window regains focus, refresh it
        case sf::Event::GainedFocus:
            param.brot->refreshWindow();
            break;

            //if the event is a resize, handle it
        case sf::Event::Resized:
            handleResize(param.brot, &param.event);
            break;

            //if the window is closed
        case sf::Event::Closed:
            param.brot->close();
            break;

        default:
            break;
    } //end event switch
}


//this function handles all keyboard input
void handleKeyboard(MandelbrotViewer *brot, sf::Event *event) {
    float color_inc = interpolate(0, 1, 30);
    switch(event->key.code) {
        //if Q, quit and close the window
        case sf::Keyboard::Q:
            brot->close();
            break;
        //if Esc, quit and close the window
        case sf::Keyboard::Escape:
            brot->close();
            break;
        //if up arrow, increase iterations
        case sf::Keyboard::Up:
            brot->incIterations();
            if (param.done) {
                handleGenerate();
                brot->updateMandelbrot();
                brot->refreshWindow();
            }
            break;
        //if down arrow, decrease iterations
        case sf::Keyboard::Down:
            brot->decIterations();
            if (param.done) {
                handleGenerate();
                brot->updateMandelbrot();
                brot->refreshWindow();
            }
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
        //if L, lock/unlock the color
        case sf::Keyboard::L:
            brot->lockColor();
            break;
        case sf::Keyboard::H:
            brot->enableOverlay(true);
            while(true) {
                brot->waitEvent(*event);
                if (event->type == sf::Event::KeyPressed) {
                    if (event->key.code == sf::Keyboard::H || event->key.code == sf::Keyboard::Escape) {
                        break;
                    } else if (event-> key.code == sf::Keyboard::Q) {
                        brot->close();
                        break;
                    }
                } else if (event->type == sf::Event::Closed) {
                    brot->close();
                }
            }
            brot->enableOverlay(false);
            break;
        //if PageUp, rotate ccw
        case sf::Keyboard::PageUp:
            brot->setRotation(brot->getRotation() + handleRotate());
            break;
        //if PageDown, rotate clockwise
        case sf::Keyboard::PageDown:
            brot->setRotation(brot->getRotation() + handleRotate());
            break;
        //if Home, set the rotation to default
        case sf::Keyboard::Home:
            brot->setRotation(0);
            break;
        //if Plus, increase oversampling level
        case sf::Keyboard::Equal:
            brot->setOversampling(brot->getOversampling()+1);
            if (param.done) {
                handleGenerate();
                brot->updateMandelbrot();
                brot->refreshWindow();
            }
            break;
        //if Minus, decrease oversampling level
        case sf::Keyboard::Dash:
            if (brot->getOversampling() > 1) 
                brot->setOversampling(brot->getOversampling()-1);
            if (param.done) {
                handleGenerate();
                brot->updateMandelbrot();
                brot->refreshWindow();
            }
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
    std::thread thread(&zoom);

    //start generating while it's zooming
    brot->generate();

    //wait for the thread to finish (wait for the zoom to finish)
    thread.join();

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
    old_center = sf::Vector2f(brot->getResWidth()/2.0, brot->getResHeight()/2.0);

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
    param.brot->setWindowActive(true);
    double inc_drag_x = interpolate(param.oldc.x, param.newc.x, param.frames);
    double inc_drag_y = interpolate(param.oldc.y, param.newc.y, param.frames);
    double inc_zoom = interpolate(1.0, param.zoom, param.frames);

    //animate the zoom
    for (int i=0; i<param.frames; i++) {
        param.newc.x = param.oldc.x + i * inc_drag_x;
        param.newc.y = param.oldc.y + i * inc_drag_y;
        param.brot->changePosView(param.newc, 1 + i * inc_zoom);
        param.brot->refreshWindow();
    }
    param.brot->setWindowActive(false);
}

//rotates the view until the key is released, then returns the new rotation
double handleRotate() {
    float rotate_inc = 1.0;
    //float rotation = param.brot->getRotation() * 180 / PI;
    float rotation = 0;
    int framerate = param.brot->getFramerate();

    //set the framerate high so that it will rotate in real time
    param.brot->setFramerate(500);

    //if page up, rotate ccw
    while (sf::Keyboard::isKeyPressed(sf::Keyboard::PageUp)) {
        param.brot->rotateView(rotation);
        param.brot->refreshWindow();
        rotation += rotate_inc;
        if (rotation >= 360) rotation -= 360;
    }

    //if page down, rotate cw
    while (sf::Keyboard::isKeyPressed(sf::Keyboard::PageDown)) {
        param.brot->rotateView(rotation);
        param.brot->refreshWindow();
        rotation -= rotate_inc;
        if (rotation < 0) rotation += 360;
    }

    param.brot->setFramerate(framerate);

    //return the rotation in radians
    return rotation * PI / 180;
}

//handle the window resize
void handleResize(MandelbrotViewer *brot, sf::Event *event) {
    int newX = event->size.width,
        newY = event->size.height;
    brot->resizeWindow(newX, newY);
    brot->generate();
    brot->updateMandelbrot();
    brot->refreshWindow();
}

void handleGenerate() {
    //start checking for events in a separate thread
    param.done = false;
    std::thread thread(eventPoll);
    thread.detach();

    //start generating
    param.brot->generate();
    param.done = true; //signal the eventPoll thread to return
}

void eventPoll() {
    while(!param.done) {
        if (param.brot->pollEvent(param.event)) {
            if (param.event.type == sf::Event::KeyPressed && 
                    (param.event.key.code == sf::Keyboard::Up    ||
                     param.event.key.code == sf::Keyboard::Down  ||
                     param.event.key.code == sf::Keyboard::Equal ||
                     param.event.key.code == sf::Keyboard::Dash)) {
                //if input is found:
                //handle the interrupt,
                handleEvent();
                //send the generation an interrupt,
                param.brot->restartGeneration();
            }
        }
        //sleep a bit so that the processor doesn't spike checking for events
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}
