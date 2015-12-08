#ifndef MANDELBROTVIEWER_H
#define MANDELBROTVIEWER_H

#define THREAD_COUNT    4

#include <SFML/Graphics.hpp>
#include <vector>

struct Color {
    int r;
    int g;
    int b;
};

class MandelbrotViewer {
    public:
        //This constructor creates a new viewer with specified resolution
        MandelbrotViewer(int resolution);
        ~MandelbrotViewer();

        //Accesor functions:
        int getResolution() {return resolution;}
        int getFramerate() {return framerateLimit;}
        double getColorMultiple() {return color_multiple;}
        sf::Vector2i getMousePosition();
        sf::Vector2f getViewCenter() {return view->getCenter();}
        sf::Vector2f getMandelbrotCenter();
        bool getEvent(sf::Event&);
        bool isOpen();
        
        //Setter functions:
        void setIterations(int iter) {max_iter = iter;}
        void setColorMultiple(double mult) {color_multiple = mult;}
        void setFramerate(int rate) {framerateLimit = rate;}
        void setColorScheme(int newScheme) {scheme = newScheme; initPalette();}
        
        //Functions to change parameters for mandelbrot generation:
        void changeColor();
        void changePos(sf::Vector2<double> new_center, double zoom_factor);
        void changePosView(sf::Vector2f new_center, double zoom_factor);

        //Functions ot generate the mandelbrot:
        void generate();

        //Functions to reset or update:
        void resetMandelbrot();
        void refreshWindow();
        void resetView();
        void close();
        void updateMandelbrot();

        //Other functions:
        void saveImage();

        //Converts a vector from pixel coordinates to the corresponding
        //coordinates of the complex plane
        sf::Vector2<double> pixelToComplex(sf::Vector2f);

    private:
        int resolution;
        int framerateLimit;
        int nextLine;

        //These are pointers to each instance's window and view
        //since we can't initialize them yet
        sf::RenderWindow *window;
        sf::View *view;

        sf::Sprite sprite;
        sf::Image image;
        sf::Texture texture;

        //Parameters to generate the mandelbrot:
        
        //this is the area of the complex plane to generate
        sf::Rect<double> area;
        
        //this changes how the colors are displayed
        double color_multiple;
        int scheme;

        //this array stores the number of iterations for each pixel
        std::vector< std::vector<int> > image_array;

        //maximum number of iterations to check for. Higher values are slower,
        //but more precise
        int max_iter;

        //Functions:
        
        //interpolate returns the increment to get from min to max in range iterations
        double interpolate(double min, double max, int range) {return (max-min)/range;}
        double interpolate(double length, int range) {return length/range;}

        //escape calculates the escape-time of given point of the mandelbrot
        int escape(double x, double y);

        //genLine is a function for worker threads: it generates the next line of the
        //mandelbrot, then moves onto the next, until the entire mandelbrot is generated
        void genLine();

        //This looks up a color to print according to the escape value given
        sf::Color findColor(int iter);

        //initialize the color palette. Having a palette helps avoid regenerating the
        //color scheme each time it is needed
        int palette[3][256];
        void initPalette();
        void smoosh(sf::Color c1, sf::Color c2, int min, int max);
};

#endif
