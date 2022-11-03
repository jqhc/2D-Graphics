#ifndef CANVAS2D_H
#define CANVAS2D_H

#include <QLabel>
#include <QMouseEvent>
#include <array>
#include "rgba.h"

class Canvas2D : public QLabel {
    Q_OBJECT
public:
    int m_width = 0;
    int m_height = 0;

    void init();
    void clearCanvas();
    bool loadImageFromFile(const QString &file);
    void displayImage();
    void resize(int w, int h);

    // This will be called when the settings have changed
    void settingsChanged();

    // Filter TODO: implement
    void filterImage();

    int posToIndex(int x, int y, int width);
    RGBA getPixelRepeated(std::vector<RGBA> &data, int x, int y);
    float getPixelRepeated(std::vector<float> &data, int x, int y);

private:
    std::vector<RGBA> m_data;

    // stores mask intensity data
    std::vector<float> mask_vector;

    // stores the colours "picked up" from canvas, for smudge brush purposes
    std::vector<RGBA> smudge_data;

    /**stores current brush type, radius and density values, to check if any
     * have been changed. if so, mask will be regenerated.*/
    int curr_brush_type;
    int curr_radius;
    int curr_density;

    void mouseDown(int x, int y);
    void mouseDragged(int x, int y);
    void mouseUp(int x, int y);

    // These are functions overriden from QWidget that we've provided
    // to prevent you from having to interact with Qt's mouse events.
    // These will pass the mouse coordinates to the above mouse functions
    // that you will have to fill in.
    virtual void mousePressEvent(QMouseEvent* event) override {
        auto [x, y] = std::array{ event->position().x(), event->position().y() };
        mouseDown(static_cast<int>(x), static_cast<int>(y));
    }
    virtual void mouseMoveEvent(QMouseEvent* event) override {
        auto [x, y] = std::array{ event->position().x(), event->position().y() };
        mouseDragged(static_cast<int>(x), static_cast<int>(y));
    }
    virtual void mouseReleaseEvent(QMouseEvent* event) override {
        auto [x, y] = std::array{ event->position().x(), event->position().y() };
        mouseUp(static_cast<int>(x), static_cast<int>(y));
    }

    // TODO: add any member variables or functions you need

    // calculates intensities and colours for each pixel on brush mask
    std::tuple<float, float, RGBA> getIntensitiesAndColour(int mask_index, int canvas_index);

    // "paints" by changing canvas colour values
    void paintBrush(int x, int y);

    // regenerates mask according to current brush type, radius, and density
    void changeMask(float (&intensity) (float));

    // checks brush type and calls changeMask with the corresponding type
    void setMask();

    // calculates the distance from a pixel on the mask to the center
    float maskDist(int x, int y);

    /**"picks up" all the colour values within the radius, at the given x and y,
     * and stores them in the vector smudge_data*/
    void storeSmudgeData(int x, int y);


    /**
     * FILTER
     */

    // different filter methods:

    // sobel edge detection
    void filterEdgeDetect(float edgeDetectSensitivity);
    // converts canvas to grayscale
    void filterGray();
    // vectors for storing intermediate canvas data
    std::vector<float> sobelXData;
    std::vector<float> sobelYData;
    std::vector<float> tempSobel;
    // methods for convolving with separated sobel kernels
    void sobelXVertical();
    void sobelXHorizontal();
    void sobelYVertical();
    void sobelYHorizontal();
    // calculates gradient and updates canvas
    void setGradient(float edgeDetectSensitivity);

    // blurring
    void filterBlur(int blurRadius);
    // convolves kernels over x and y directions
    void convolveX(const std::vector<float> &kernel, int radius);
    void convolveY(const std::vector<float> &kernel, int radius);
    // vector for storing blur kernel (triangle)
    std::vector<float> blur_kernel_1D;
    // updates blur kernel for given radius
    void setBlurKernel(int blurRadius);

    // scaling
    void filterScale(float scaleX, float scaleY);
    // vector for storing intermediate canvas data
    std::vector<RGBA> output;
    // scales in x and y directions
    void scaleXDirection(float a);
    void scaleYDirection(float a);
    // backmapping function for getting input coordinate from output coordinate
    float backMap(int x, float &a);
    // calculates weights of triangle filter
    float triangleWeight(float distance, float &radius);
    // performs the resampling convolution over x and y directions
    RGBA h_prime_x(int x, int y, float &a, float &radius);
    RGBA h_prime_y(int x, int y, float &a, float &radius);

    // median
    void filterMedian(int medianRadius);
    RGBA calculateMedian(int x, int y, int medianRadius);
};

#endif // CANVAS2D_H
