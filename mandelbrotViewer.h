#ifndef MANDELBROTVIEWER_H
#define MANDELBROTVIEWER_H

#include <SFML/Graphics.hpp>
#include <gmp.h>
#include <vector>

struct Color {
    int r;
    int g;
    int b;
};

class MandelbrotViewer {
    public:
        //This constructor creates a new viewer with specified resolution
        MandelbrotViewer(int res_x, int res_y);
        ~MandelbrotViewer();

        //Accesor functions:
        int getResWidth() {return res_width;}
        int getResHeight() {return res_height;}
        int getFramerate() {return framerateLimit;}
        int getIters() {return max_iter;}
        double getRotation() {return rotation;}
        double getColorMultiple() {return color_multiple;}
        sf::Vector2i getMousePosition();
        sf::Vector2f getViewCenter() {return view->getCenter();}
        sf::Vector2f getMandelbrotCenter();
        bool waitEvent(sf::Event&);
        bool pollEvent(sf::Event&);
        bool isColorLocked() {return color_locked;}
        bool isOpen();
        
        //Setter functions:
        void incIterations();
        void decIterations();
        void setIterations(int iter) {max_iter = iter; initPalette();}
        void setColorMultiple(double mult) {color_multiple = mult;}
        void setFramerate(int rate) {framerateLimit = rate;}
        void setColorScheme(int newScheme);
        void setRotation(double radians);
        void restartGeneration() {restart_gen = true;}
        void lockColor();
        
        //Functions to change parameters for mandelbrot generation:
        void changeColor();
        void changePos(sf::Vector2<double> new_center, double zoom_factor);
        void changePosView(sf::Vector2f new_center, double zoom_factor);
        void resizeWindow(int newX, int newY);

        //Functions to generate the mandelbrot:
        void generate();

        //Functions to reset or update:
        void resetMandelbrot();
        void refreshWindow();
        void resetView();
        void close();
        void updateMandelbrot();
        void setWindowActive(bool);

        //Other functions:
        void saveImage(); //save the image as a png in the local folder
        void enableOverlay(bool); //enable a help overlay with controls, etc.
        void rotateView(float angle);

        //Converts a vector from pixel coordinates to the corresponding
        //coordinates of the complex plane
        sf::Vector2<double> pixelToComplex(sf::Vector2f);

    private:
        int res_height;
        int res_width;
        int framerateLimit;
        int nextLine;

        sf::Sprite sprite;
        sf::Image image;
        sf::Texture texture;
        sf::Font font;

        //These are pointers to each instance's window and view
        //since we can't initialize them yet
        sf::RenderWindow *window;
        sf::View *view;

        //Parameters to generate the mandelbrot:
        bool restart_gen; //set to true to stop generation before it's finished

        //this is the area of the complex plane to generate
        sf::Rect<mpf_t> area;
        mpf_t area_inc; //this is complex plane area per pixel

        //this is the current rotation of the mandelbrot - 0 radians is positive x axis
        double rotation;
        
        //this changes how the colors are displayed
        double color_multiple;
        bool color_locked;
        int scheme;

        //Holds the maximum number of concurrent threads suppported by the current CPU
        unsigned int max_threads;

        //this array stores the number of iterations for each pixel
        std::vector< std::vector<int> > image_array;

        //maximum number of iterations to check for. Higher values are slower,
        //but more precise
        unsigned int max_iter;
        unsigned int last_max_iter;


        //Functions:
        
        //interpolate returns the increment to get from min to max in range iterations
        double interpolate(double min, double max, int range) {return (max-min)/range;}
        double interpolate(double length, int range) {return length/range;}

        //escape calculates the escape-time of given point of the mandelbrot
        int escape(int row, int column);

        //genLine is a function for worker threads: it generates the next line of the
        //mandelbrot, then moves onto the next, until the entire mandelbrot is generated
        void genLine();

        //this looks up a color to print according to the escape value given
        sf::Color findColor(int iter);

        //this function handles rotation - it takes in a complex point with zero rotation
        //and returns where that point is when rotated
        sf::Vector2<double> rotate(sf::Vector2<double>);

        //initialize the color palette. Having a palette helps avoid regenerating the
        //color scheme each time it is needed
        std::vector< std::vector<int> > palette;
        void initPalette();
        void smoosh(sf::Color c1, sf::Color c2, float min, float max);
};

#endif
