#include "mandelbrotViewer.h"
#include <iostream>


//this struct holds the parameters for the zoom function,
//since it needs to be threaded and cannot accept arguments
struct zoomParameters {
    MandelbrotViewer *viewer;
    sf::Vector2f oldc;
    sf::Vector2f newc;
    double zoom;
    int frames;
} param;

//returns range increments to get from min to max
double interpolate(double min, double max, int range) { return (max-min)/range; }

//animates the zoom of the viewer (so the viewer zooms while the mandelbrot is generating)
//uses the struct param as parameters
void zoom() {
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
}

int main(int argc, char **argv) {
    if (argc > 2) {
        std::cout << "Doing a fixed test - remember to time!" << std::endl;
        std::cout << "Also remember to call like: " << argv[0] << " <#iterations> <zoom factor>" << std::endl;
        MandelbrotViewer brot(1024);
        brot.resetMandelbrot();
        brot.setIterations(atoi(argv[1]));
        sf::Vector2<double> new_pos;
        new_pos.x = 0.013438870532012129028364919004019686867528573314565492885548699;
        new_pos.y = 0.655614218769465062251320027664617466691295975864786403994151735;
        double zoom = atof(argv[2]);
        brot.changePos(new_pos, zoom);
        brot.generate();
        brot.updateMandelbrot();
        brot.saveImage();
        return 0;
    }

    int resolution = 720;
    int iterations = 100;
    int framerateLimit;
    float color_inc = interpolate(0, 1, 30);
    
    //get vectors ready
    sf::Vector2f old_position;
    sf::Vector2f new_position;
    sf::Vector2f old_center;
    sf::Vector2f new_center;
    sf::Vector2f difference;
    sf::Vector2i temp;

    //create the mandelbrotviewer instance
    MandelbrotViewer brot(resolution);
    brot.resetMandelbrot();
    brot.generate();
    brot.updateMandelbrot();
    brot.refreshWindow();

    sf::Event event;
    framerateLimit = brot.getFramerate();

    //point the zoom function to the 'brot' instance
    param.viewer = &brot;

    //main window loop
    while (brot.isOpen()) {

        //main event loop
        while (brot.getEvent(event)) {

            //this big switch statement handles all types of input
            switch (event.type) {

                //if the event is a keypress
                case sf::Event::KeyPressed:

                    //another switch statement handles different keys
                    switch (event.key.code) {

                        //if Q, quit and close the window
                        case sf::Keyboard::Q:
                            brot.close();
                            break;
                        //if up arrow, increase iterations
                        case sf::Keyboard::Up:
                            iterations += 30;
                            brot.setIterations(iterations);
                            brot.generate();
                            brot.updateMandelbrot();
                            brot.refreshWindow();
                            break;
                        //if down arrow, decrease iterations
                        case sf::Keyboard::Down:
                            iterations -= 30;
                            if (iterations < 100) iterations = 100;
                            brot.setIterations(iterations);
                            brot.generate();
                            brot.updateMandelbrot();
                            brot.refreshWindow();
                            break;
                        //if right arrow, increase color_multiple until released
                        case sf::Keyboard::Right:
                            color_inc = interpolate(0, 1, 25);
                            while (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) {
                                brot.setColorMultiple(brot.getColorMultiple() + color_inc);
                                brot.changeColor();
                                brot.updateMandelbrot();
                                brot.refreshWindow();
                            }
                            break;
                        //if left arrow, decrease color_multiple until released
                        case sf::Keyboard::Left:
                            color_inc = interpolate(1, 0, 25);
                            while (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) {
                                if (brot.getColorMultiple() > 1) {
                                    brot.setColorMultiple(brot.getColorMultiple() + color_inc);
                                    brot.changeColor();
                                    brot.updateMandelbrot();
                                    brot.refreshWindow();
                                }
                            }
                            break;
                        //if it's a 1, change to color scheme 1
                        case sf::Keyboard::Num1:
                            brot.setColorScheme(1);
                            brot.changeColor();
                            brot.updateMandelbrot();
                            brot.refreshWindow();
                            break;
                        //if it's a 2, change to color scheme 2
                        case sf::Keyboard::Num2:
                            brot.setColorScheme(2);
                            brot.changeColor();
                            brot.updateMandelbrot();
                            brot.refreshWindow();
                            break;
                        //if it's a 3, change to color scheme 3
                        case sf::Keyboard::Num3:
                            brot.setColorScheme(3);
                            brot.changeColor();
                            brot.updateMandelbrot();
                            brot.refreshWindow();
                            break;
                        //if R, reset the mandelbrot to the starting image
                        case sf::Keyboard::R:
                            brot.resetMandelbrot();
                            brot.resetView();
                            brot.generate();
                            brot.updateMandelbrot();
                            brot.refreshWindow();
                            break;
                        //if S, save the current image
                        case sf::Keyboard::S:
                            brot.saveImage();
                            break;
		    case sf::Keyboard::M:
	    case sf::Keyboard::N:
		  {
                    //set the zoom frames up so that it is more animated
                    param.frames = 30;
                    old_center = brot.getViewCenter();
                    new_center.x = old_center.x;//event.mouseWheelScroll.x;
                    new_center.y = old_center.y;//event.mouseWheelScroll.y;

                    //if it's an upward scroll, get ready to zoom in
                    //if (event.mouseWheelScroll.delta > 0) {
		    if (event.key.code == sf::Keyboard::M) {
                        brot.changePos(brot.pixelToComplex(new_center), 0.5);
                        param.zoom = 0.5;
                    } //if it's a downward scroll, get ready to zoom out
                    else {
                        brot.changePos(brot.pixelToComplex(new_center), 2.0);
                        param.zoom = 2.0;
                    }

                    //set the parameters for the zoom
                    param.oldc = old_center;
                    param.newc = new_center;

                    //start zooming with a worker thread, so that it can generate
                    //the new image while it's zooming
                    sf::Thread thread(&zoom);
                    thread.launch();
                    
                    //start generating while it's zooming
                    brot.generate();

                    //wait for the thread to finish (wait for the zoom to finish)
                    thread.wait();

                    //now display the new mandelbrot
                    brot.updateMandelbrot();
                    brot.resetView();
                    brot.refreshWindow();
                    break; 
                }
                        //end the keypress switch
                        default:
                            break;
                    }
                    //end the keypress case
                    break;

                //if the event is a mousewheel scroll, zoom in with the
                //mouse coordinates as the new center
                //NOTE: this has to have brackets because it is declaring new variables
#if 0
		    //                case sf::Event::MouseWheelScrolled:
	    case sf::Event::Keyboard::M:
	    case sf::Event::Keyboard::N:
		  {
                    //set the zoom frames up so that it is more animated
                    param.frames = 30;
                    old_center = brot.getViewCenter();
                    //new_center.x = event.mouseWheelScroll.x;
                    //new_center.y = event.mouseWheelScroll.y;

                    //if it's an upward scroll, get ready to zoom in
                    if (event.mouseWheelScroll.delta > 0) {
                        brot.changePos(brot.pixelToComplex(new_center), 0.5);
                        param.zoom = 0.5;
                    } //if it's a downward scroll, get ready to zoom out
                    else if (event.mouseWheelScroll.delta < 0) {
                        brot.changePos(brot.pixelToComplex(new_center), 2.0);
                        param.zoom = 2.0;
                    }

                    //set the parameters for the zoom
                    param.oldc = old_center;
                    param.newc = new_center;

                    //start zooming with a worker thread, so that it can generate
                    //the new image while it's zooming
                    sf::Thread thread(&zoom);
                    thread.launch();
                    
                    //start generating while it's zooming
                    brot.generate();

                    //wait for the thread to finish (wait for the zoom to finish)
                    thread.wait();

                    //now display the new mandelbrot
                    brot.updateMandelbrot();
                    brot.resetView();
                    brot.refreshWindow();
                    break; 
                }
#endif

                //if the event is a click, drag the view:
                case sf::Event::MouseButtonPressed:

                    //set drag frames low so that it will move in real time
                    param.frames = 2;

                    //save the old center/mouse_position as reference
                    temp = brot.getMousePosition();
                    old_position.x = temp.x;
                    old_position.y = temp.y;
                    old_center = brot.getViewCenter();
                    param.oldc.x = old_center.x;
                    param.oldc.y = old_center.y;

                    //set zoom to 1 so that it doesn't change, only drags
                    param.zoom = 1.0;

                    //set the framrate very high so that it will drag in real time
                    brot.setFramerate(500);

                    while (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {

                        //get the new mouse position
                        temp = brot.getMousePosition();
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
                    brot.setFramerate(framerateLimit);

                    //now regenerate the mandelbrot at the new position.
                    //This can't be threaded like zoom, because it doesn't
                    //know where to generate at until it's done dragging
                    temp = brot.getMousePosition();
                    new_position.x = temp.x;
                    new_position.y = temp.y;

                    //get the difference of how far it's moved
                    difference.x = new_position.x - old_position.x;
                    difference.y = new_position.y - old_position.y;

                    //set the old center
                    old_center = sf::Vector2f(resolution/2, resolution/2);

                    //calculate the new center
                    new_center = old_center - difference;

                    brot.changePos(brot.pixelToComplex(new_center), 1.0);
                    brot.generate();
                    brot.resetView();
                    brot.updateMandelbrot();
                    brot.refreshWindow();
                    break;
                
                default:
                    break;

            } //end the main event switch statement

        } //end the main event loop

    } //end main window loop

    return 0;
} //end main
