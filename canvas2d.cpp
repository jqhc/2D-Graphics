#include "canvas2d.h"
#include <QPainter>
#include <QMessageBox>
#include <QFileDialog>
#include <iostream>
#include "settings.h"
#include <cmath>
#include <stdlib.h>
#include <list>

/**
 * @brief Initializes new 500x500 canvas
 */
void Canvas2D::init() {
    m_width = 500;
    m_height = 500;
    clearCanvas();

    // set initial brush type to linear
    curr_brush_type = settings.brushType;
    curr_radius = settings.brushRadius;
    curr_density = settings.brushDensity;

    // initiate mask
    setMask();

    // reserve space in sobel data vectors
    tempSobel.reserve(m_width*m_height);
    sobelXData.reserve(m_width*m_height);
    sobelYData.reserve(m_width*m_height);
}

/**
 * @brief Canvas2D::clearCanvas sets all canvas pixels to blank white
 */
void Canvas2D::clearCanvas() {
    m_data.assign(m_width * m_height, RGBA{255, 255, 255, 255});
    settings.imagePath = "";
    displayImage();
}

/**
 * @brief Stores the image specified from the input file in this class's
 * `std::vector<RGBA> m_image`.
 * Also saves the image width and height to canvas width and height respectively.
 * @param file: file path to an image
 * @return True if successfully loads image, False otherwise.
 */
bool Canvas2D::loadImageFromFile(const QString &file) {
    QImage myImage;
    if (!myImage.load(file)) {
        std::cout<<"Failed to load in image"<<std::endl;
        return false;
    }
    myImage = myImage.convertToFormat(QImage::Format_RGBX8888);
    m_width = myImage.width();
    m_height = myImage.height();
    QByteArray arr = QByteArray::fromRawData((const char*) myImage.bits(), myImage.sizeInBytes());

    m_data.clear();
    m_data.reserve(m_width * m_height);
    for (int i = 0; i < arr.size() / 4.f; i++){
        m_data.push_back(RGBA{(std::uint8_t) arr[4*i], (std::uint8_t) arr[4*i+1], (std::uint8_t) arr[4*i+2], (std::uint8_t) arr[4*i+3]});
    }

    displayImage();
    return true;
}

/**
 * @brief Get Canvas2D's image data and display this to the GUI
 */
void Canvas2D::displayImage() {
    QByteArray* img = new QByteArray(reinterpret_cast<const char*>(m_data.data()), 4*m_data.size());
    QImage now = QImage((const uchar*)img->data(), m_width, m_height, QImage::Format_RGBX8888);
    setPixmap(QPixmap::fromImage(now));
    setFixedSize(m_width, m_height);
    update();
}

/**
 * @brief Canvas2D::resize resizes canvas to new width and height
 * @param w
 * @param h
 */
void Canvas2D::resize(int w, int h) {
    m_width = w;
    m_height = h;
    m_data.resize(w * h);
    displayImage();
}

/**
 * @brief Called when the filter button is pressed in the UI
 */
void Canvas2D::filterImage() {
    switch(settings.filterType) {
        case FILTER_EDGE_DETECT:
            filterEdgeDetect(settings.edgeDetectSensitivity);
            break;
        case FILTER_BLUR:
            filterBlur(settings.blurRadius);
            break;
        case FILTER_SCALE:
            filterScale(settings.scaleX, settings.scaleY);
            break;
        case FILTER_MEDIAN:
            filterMedian(settings.medianRadius);
            break;
        default:
            break;
    }
    displayImage();
}

void Canvas2D::filterEdgeDetect(float edgeDetectSensitivity) {
    // first convert image to grayscale
    filterGray();

    // reserve space in temporary, G_x and G_y vectors
    tempSobel.reserve(m_width*m_height);
    sobelXData.reserve(m_width*m_height);
    sobelYData.reserve(m_width*m_height);

    // populate G_x
    sobelXHorizontal();
    sobelXVertical();

    // populate G_y
    sobelYHorizontal();
    sobelYVertical();

    // calculate G and set canvas data
    setGradient(edgeDetectSensitivity);
}

std::uint8_t rgbaToGray(const RGBA &pixel) {
    // multiplies each colour value by a coefficient and adds them to get a
    // grayscale value
    return (0.299*pixel.r + 0.587*pixel.g + 0.114*pixel.b);
}

void Canvas2D::filterGray() {
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            RGBA &currentPixel = m_data[posToIndex(x, y, m_width)];

            // convert pixel to grayscale value
            std::uint8_t grayPixel = rgbaToGray(currentPixel);
            // update pixel values
            currentPixel.r = grayPixel;
            currentPixel.g = grayPixel;
            currentPixel.b = grayPixel;
        }
    }
}

void Canvas2D::sobelXHorizontal() {
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            // calculate new intensity with weights:
            // -1 0 1
            float newIntensity =
                    (float)getPixelRepeated(m_data, x+1, y).r // right pixel
                    - (float)getPixelRepeated(m_data, x-1, y).r; // left pixel
            // store in temp canvas
            tempSobel[posToIndex(x,y,m_width)] = newIntensity;
        }
    }
    // copy temp pixels to sobel x vector
    for (int i = 0; i < m_width*m_height; i++) {
        sobelYData[i] = tempSobel[i];
    }
}

void Canvas2D::sobelXVertical() {
    for (int x = 0; x < m_width; x++) {
        for (int y = 0; y < m_height; y++) {
            // calculate new intensity with weights:
            // 1
            // 2
            // 1
            float newIntensity =
                    getPixelRepeated(sobelXData, x, y-1) // up pixel
                    + 2.0f*getPixelRepeated(sobelXData, x, y) // current pixel
                    + getPixelRepeated(sobelXData, x, y+1); // down pixel
            // store in temp canvas
            tempSobel[posToIndex(x,y,m_width)] = newIntensity;
        }
    }
    // copy temp pixels to sobel x vector
    for (int i = 0; i < m_width*m_height; i++) {
        sobelYData[i] = tempSobel[i];
    }
}

void Canvas2D::sobelYHorizontal() {
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            // calculate new intensity with weights:
            // 1 2 1
            float newIntensity =
                    (float)getPixelRepeated(m_data, x-1, y).r // up pixel
                    + 2.0f*(float)getPixelRepeated(m_data, x, y).r // current pixel
                    + (float)getPixelRepeated(m_data, x+1, y).r; // down pixel
            // store in temp canvas
            tempSobel[posToIndex(x,y,m_width)] = newIntensity;
        }
    }
    // copy temp pixels to sobel y vector
    for (int i = 0; i < m_width*m_height; i++) {
        sobelYData[i] = tempSobel[i];
    }
}

void Canvas2D::sobelYVertical() {
    for (int x = 0; x < m_width; x++) {
        for (int y = 0; y < m_height; y++) {
            // calculate new intensity with weights:
            // -1
            // 0
            // 1
            float newIntensity =
                    getPixelRepeated(sobelYData, x, y+1) // down pixel
                    - getPixelRepeated(sobelYData, x, y-1); // up pixel
            // store in temp canvas
            tempSobel[posToIndex(x,y,m_width)] = newIntensity;
        }
    }
    // copy temp pixels to sobel y vector
    for (int i = 0; i < m_width*m_height; i++) {
        sobelYData[i] = tempSobel[i];
    }
}


void Canvas2D::setGradient(float edgeDetectSensitivity) {

    for (int i = 0; i < m_height*m_width; i++) {
        // get G_x and G_y values
        float g_x_pixel = sobelXData[i];
        float g_y_pixel = sobelYData[i];

        // calculate magnitude of gradient
        float g_pixel = sqrt(pow(g_x_pixel,2) + pow(g_y_pixel, 2));

        // multiply by sensitivity value
        g_pixel *= edgeDetectSensitivity;

        // clamp to uint8_t
        std::uint8_t rgba_g = (g_pixel <= 255) ? floor(g_pixel+0.5f) : 255;

        // set canvas pixel to this value
        m_data[i] = RGBA{rgba_g, rgba_g, rgba_g, 255};
    }
}

void Canvas2D::filterBlur(int blurRadius) {
    setBlurKernel(blurRadius);
    convolveX(blur_kernel_1D, blurRadius);
    convolveY(blur_kernel_1D, blurRadius);
}

void Canvas2D::setBlurKernel(int blurRadius) {
    int kernel_diameter = 2 * blurRadius + 1;

    // clear current kernel
    blur_kernel_1D.clear();
    // reserve 2R+1 elements in the kernel vector
    blur_kernel_1D.reserve(kernel_diameter);

    // populate the triangle kernel vector
    int r;
    for (r = 1; r <= blurRadius + 1; r++) {
        blur_kernel_1D.push_back(r);
    }
    for (r = blurRadius; r >0; r--) {
        blur_kernel_1D.push_back(r);
    }

    // define the factor that normalizes the sum of the vector to 1.
    float normalizing_factor = pow(blurRadius+1, 2);

    for (int r = 0; r < kernel_diameter; r++) {
        blur_kernel_1D[r] /= normalizing_factor;
    }
}

/**
 * @brief Given a 1D kernel and its radius, convolves it horizontally over each row
 *        of the image and updates it.
 */
void Canvas2D::convolveX(const std::vector<float> &kernel, int radius) {

    // temporary vector to store convolved pixels
    std::vector<RGBA> result(m_width * m_height);

    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            // set up accumulators for each colour
            float redAcc = 0;
            float greenAcc = 0;
            float blueAcc = 0;

            for (int kern_x = 0; kern_x < 2*radius+1; kern_x++) {
                // Get kernel weight
                const float &weight = kernel[kern_x];
                // Get canvas pixel
                const RGBA &canvas_pixel = getPixelRepeated(m_data, x - radius + kern_x, y);

                // add multiplied values to accumulators
                redAcc += weight * (float)canvas_pixel.r;
                greenAcc += weight * (float)canvas_pixel.g;
                blueAcc += weight * (float)canvas_pixel.b;
            }

            // set and create new pixel
            std::uint8_t newRed = floor(redAcc + 0.5f);
            std::uint8_t newGreen = floor(greenAcc + 0.5f);
            std::uint8_t newBlue = floor(blueAcc + 0.5f);
            result[posToIndex(x,y,m_width)] = RGBA{newRed, newGreen, newBlue, 255};
        }
    }
    // copy result pixels to canvas
    for (int i = 0; i < m_width*m_height; i++) {
        m_data[i] = result[i];
    }
}

/**
 * @brief Given a 1D kernel and its radius, convolves it vertically over each column
 *        of the image and updates it.
 */
void Canvas2D::convolveY(const std::vector<float> &kernel, int radius) {

    // temporary vector to store convolved pixels
    std::vector<RGBA> result(m_width * m_height);

    for (int x = 0; x < m_width; x++) {
        for (int y = 0; y < m_height; y++) {
            float redAcc = 0;
            float greenAcc = 0;
            float blueAcc = 0;

            for (int kern_y = 0; kern_y < 2*radius+1; kern_y++) {
                // Get kernel weight
                const float &weight = kernel[kern_y];
                // Get canvas pixel
                const RGBA &canvas_pixel = getPixelRepeated(m_data, x, y - radius + kern_y);

                redAcc += weight * (float)canvas_pixel.r;
                greenAcc += weight * (float)canvas_pixel.g;
                blueAcc += weight * (float)canvas_pixel.b;
            }

            std::uint8_t newRed = floor(redAcc + 0.5f);
            std::uint8_t newGreen = floor(greenAcc + 0.5f);
            std::uint8_t newBlue = floor(blueAcc + 0.5f);
            result[posToIndex(x,y,m_width)] = RGBA{newRed, newGreen, newBlue, 255};
        }
    }
    // copy result pixels to canvas
    for (int i = 0; i < m_width*m_height; i++) {
        m_data[i] = result[i];
    }
}

/**
 * @brief scales the image in the X direction by a factor of scaleX, then scales
 *        in the Y direction by a factor of scaleY.
 */
void Canvas2D::filterScale(float scaleX, float scaleY) {
    if (scaleX != 1) {
        scaleXDirection(scaleX);
    }
    if (scaleY != 1) {
        scaleYDirection(scaleY);
    }
}

/**
 * @brief scales the image in the X direction by a factor of a.
 */
void Canvas2D::scaleXDirection(float a) {
    // calculate width of output image
    int new_width = a*m_width;

    // prepare vector for output image
    output.clear();
    output.reserve(new_width*m_height);

    // calculate radius of triangle filter
    float radius = (a > 1) ? 1 : 1/a;

    // recurse through output image, get its value by resampling from input image
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < new_width; x++) {
            // get new colour value by resampling
            RGBA newColour = h_prime_x(x, y, a, radius);

            // set result canvas pixel to new colour
            output[posToIndex(x,y,new_width)] = newColour;
        }
    }

    // change canvas width and copy output pixels to canvas
    m_width = new_width;
    m_data.resize(m_width*m_height);
    for (int i = 0; i < m_width*m_height; i++) {
        m_data[i] = output[i];
    }
}

/**
 * @brief scales the image in the Y direction by a factor of a.
 */
void Canvas2D::scaleYDirection(float a) {
    // calculate width of output image
    int new_height = a*m_height;

    // prepare vector for output image
    output.clear();
    output.reserve(m_width*new_height);

    // calculate radius of triangle filter
    float radius = (a > 1) ? 1 : 1/a;

    // recurse through output image, get its value by resampling from input image
    for (int x = 0; x < m_width; x++) {
        for (int y = 0; y < new_height; y++) {
            // get new colour value by resampling
            RGBA newColour = h_prime_y(x, y, a, radius);

            // set result canvas pixel to new colour
            output[posToIndex(x,y,m_width)] = newColour;
        }
    }

    // change canvas height and copy output pixels to canvas
    m_height = new_height;
    m_data.resize(m_width*m_height);
    for (int i = 0; i < m_width*m_height; i++) {
        m_data[i] = output[i];
    }
}

/**
 * @brief given a coordinate x in the output image, and a scaling factor a, computes
 *        the position of x in the input image.
 */
float Canvas2D::backMap(int x, float &a) {
    return ((float)x/a) + (1-a)/(2*a);
}

/**
 * @brief calculates the weight of a pixel of given distance from the center of the
 *        triangle filter with given radius.
 */
float Canvas2D::triangleWeight(float distance, float &radius) {
    // radius is either 1 or 1/a, given by caller.
    if (distance < -radius || distance > radius) {
        return 0;
    } else {
        return (1 - std::abs(distance)/radius)/radius;
    }
}

/**
 * @brief calculates canvas RGBA value at position x in output image, where x is scaled.
 * @param x - the x-coordinate of the pixel in the output image.
 * @param y - the y-coordinate of the pixel in the output image.
 * @param a - the scaling factor.
 * @param radius - the radius of the kernel.
 */
RGBA Canvas2D::h_prime_x(int x, int y, float &a, float &radius) {
    // initialize colour and weightSum accumulators
    float redAcc = 0;
    float greenAcc = 0;
    float blueAcc = 0;

    float weightSum = 0;

    // find position of x in input image
    float center = backMap(x, a);

    // get leftmost and rightmost input indices in kernel range
    int left = ceil(center - radius);
    int right = floor(center + radius);

    for (int index = left; index <= right; index++) {
        // get pixel at this index
        RGBA input_pixel = getPixelRepeated(m_data, index, y);

        // calculate weight according to triangle filter
        float weight = triangleWeight(index - center, radius);

        // add to weight sum and accumulators
        weightSum += weight;

        redAcc += weight * (float)input_pixel.r;
        greenAcc += weight * (float)input_pixel.g;
        blueAcc += weight * (float)input_pixel.b;
    }

    // normalize and cast colours with weight_sum
    std::uint8_t newRed = floor(redAcc / weightSum + 0.5f);
    std::uint8_t newGreen = floor(greenAcc / weightSum + 0.5f);
    std::uint8_t newBlue = floor(blueAcc / weightSum + 0.5f);

    // return the new colour, to be stored in the output image.
    return RGBA{newRed, newGreen, newBlue, 255};
}

/**
 * @brief calculates canvas RGBA value at position y in output image, where y is scaled.
 * @param x - the x-coordinate of the pixel in the output image.
 * @param y - the y-coordinate of the pixel in the output image.
 * @param a - the scaling factor.
 * @param radius - the radius of the kernel.
 */
RGBA Canvas2D::h_prime_y(int x, int y, float &a, float &radius) {
    // initialize colour and weightSum accumulators
    float redAcc = 0;
    float greenAcc = 0;
    float blueAcc = 0;

    float weightSum = 0;

    // find position of x in input image
    float center = backMap(y, a);

    // get leftmost and rightmost input indices in kernel range
    int top = ceil(center - radius);
    int bottom = floor(center + radius);

    for (int index = top; index <= bottom; index++) {
        // get pixel at this index
        RGBA input_pixel = getPixelRepeated(m_data, x, index);

        // calculate weight according to triangle filter
        float weight = triangleWeight(index - center, radius);

        // add to weight sum and accumulators
        weightSum += weight;

        redAcc += weight * (float)input_pixel.r;
        greenAcc += weight * (float)input_pixel.g;
        blueAcc += weight * (float)input_pixel.b;
    }

    // normalize and cast colours with weight_sum
    std::uint8_t newRed = floor(redAcc / weightSum + 0.5f);
    std::uint8_t newGreen = floor(greenAcc / weightSum + 0.5f);
    std::uint8_t newBlue = floor(blueAcc / weightSum + 0.5f);

    // return the new colour, to be stored in the output image.
    return RGBA{newRed, newGreen, newBlue, 255};
}

void Canvas2D::filterMedian(int medianRadius) {  
    // initialize temporary result vector
    std::vector<RGBA> result;
    result.reserve(m_width*m_height);

    // store median values in result
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            result[posToIndex(x, y, m_width)] = calculateMedian(x, y, medianRadius);
        }
    }

    // copy result pixels over to canvas
    for (int i = 0; i < m_width*m_height; i++) {
        m_data[i] = result[i];
    }
}

/**
 * @brief given a vector sorted from least to greatest and a new value,
 *        inserts the value while preserving sorting.
 */
void insertSorted(std::vector<uint8_t> &vals, uint8_t newColour) {
    for (std::vector<uint8_t>::iterator it = vals.begin(); it != vals.end(); it++) {
        if (*it >= newColour) {
            vals.insert(it, newColour);
            return;
        }
    }
    vals.push_back(newColour);
}

RGBA Canvas2D::calculateMedian(int x, int y, int medianRadius) {
    // calculate diameter of "kernel"
    int medianDiameter = 2*medianRadius + 1;

    // initialize vectors for storing red, green, and blue values
    std::vector<uint8_t> redVals;
    std::vector<uint8_t> greenVals;
    std::vector<uint8_t> blueVals;
    redVals.reserve(pow(medianDiameter, 2));
    greenVals.reserve(pow(medianDiameter, 2));
    blueVals.reserve(pow(medianDiameter, 2));

    // recurses through pixels within kernel and adds the colour values to the vectors
    for (int kern_y = 0; kern_y < medianDiameter; kern_y++) {
        for (int kern_x = 0; kern_x < medianDiameter; kern_x++) {
            RGBA canvas_pixel = getPixelRepeated(m_data, x - medianRadius + kern_x, y - medianRadius + kern_y);
            insertSorted(redVals, canvas_pixel.r);
            insertSorted(greenVals, canvas_pixel.g);
            insertSorted(blueVals, canvas_pixel.b);
        }
    }

    // selects median from each vector
    uint8_t medianRed = redVals[medianRadius*(medianDiameter + 1)];
    uint8_t medianGreen = greenVals[medianRadius*(medianDiameter + 1)];
    uint8_t medianBlue = blueVals[medianRadius*(medianDiameter + 1)];

    return RGBA{medianRed, medianGreen, medianBlue, 255};
}

/**
 * @brief Called when any of the parameters in the UI are modified.
 */
void Canvas2D::settingsChanged() {
    // this saves your UI settings locally to load next time you run the program
    settings.saveSettings();

    // only regenerate mask if the brush type or radius has been changed
    if (settings.brushType != curr_brush_type || settings.brushRadius != curr_radius) {
        // change mask to new brushType/radius
        setMask();
        // change current brushType and radius values
        curr_brush_type = settings.brushType;
        curr_radius = settings.brushRadius;
    // also regenerate if brush is spray and density gets changed
    } else if (settings.brushType == 4 && settings.brushDensity != curr_density) {
        // regenerate mask for new density
        setMask();
        // change current density
        curr_density = settings.brushDensity;
    }
}

/**
 * @brief Converts coordinates of an image to an index in the canvas array.
 */
int Canvas2D::posToIndex(int x, int y, int width) {
    return width* y + x;
}
/**
 * @brief Repeats the pixel on the edge of the image such that A,B,C,D looks like ...A,A,A,B,C,D,D,D...
 */
RGBA Canvas2D::getPixelRepeated(std::vector<RGBA> &data, int x, int y) {
    int newX = (x < 0) ? 0 : std::min(x, m_width  - 1);
    int newY = (y < 0) ? 0 : std::min(y, m_height - 1);
    return data[m_width * newY + newX];
}
float Canvas2D::getPixelRepeated(std::vector<float> &data, int x, int y) {
    int newX = (x < 0) ? 0 : std::min(x, m_width  - 1);
    int newY = (y < 0) ? 0 : std::min(y, m_height - 1);
    return data[m_width * newY + newX];
}

/**
 * @brief Given a current set of coordinates within the mask, returns the distance
 *        from these coordinates to the center of the mask.
 */
float Canvas2D::maskDist(int x, int y) {
    return sqrt(pow(x-settings.brushRadius, 2) + pow(y-settings.brushRadius, 2));
}

/**
 * @brief returns a tuple containing the intensity of the mask and canvas, and the colour, at the given indices.
 */
std::tuple<float, float, RGBA> Canvas2D::getIntensitiesAndColour(int mask_index, int canvas_index) {
    float mask_intensity;
    float canvas_intensity;
    RGBA curr_colour;

    if (settings.brushType == 3) { // if smudge brush, use colours from stored smudge data

        // find the intensity of the mask at this point - do not consider alpha
        mask_intensity = mask_vector[mask_index];
        // calculate the opacity of the canvas to then be 1 - mask_intensity
        canvas_intensity = 1 - mask_intensity;

        // colour to be painted is the colour from the smudge data
        curr_colour = smudge_data[mask_index];

    } else { // if not smudge brush, just use brush colour

        // find the intensity of the mask at this point
        mask_intensity = ((float)settings.brushColor.a / (float)255) * mask_vector[mask_index];

        if (settings.brushType == 4) { // if spray brush, randomly multiply intensity by 0 or 1

            // generates a random integer in the interval [1,12000] and returns it.
            int rand_int = arc4random_uniform(400) + 1;
            // multiplies by 1 if the random integer is less than the density, otherwise 0.
            mask_intensity = mask_intensity*(rand_int <= settings.brushDensity);
        }

        // calculate the opacity of the canvas to then be 1 - mask_intensity
        canvas_intensity = 1 - mask_intensity;

        // colour to be painted is brush colour
        curr_colour = settings.brushColor;
    }

    return std::tuple{mask_intensity, canvas_intensity, curr_colour};
}

/**
 * @brief When called, paints on the canvas according to the current brush settings.
 */
void Canvas2D::paintBrush(int x, int y) {
    int r = settings.brushRadius;
    int mask_diameter = 2 * r + 1;

    for (int y_iter = 0; y_iter < mask_diameter; y_iter++) {
        for (int x_iter = 0; x_iter < mask_diameter; x_iter++) {

            // finds the coordinates of (x_iter, y_iter) on the canvas
            int canvas_x = x_iter + x - r;
            int canvas_y = y_iter + y - r;

            // only paints if these canvas coordinates are within the canvas
            if (canvas_x >= 0 && canvas_x < m_width && canvas_y >= 0 && canvas_y < m_height) {

                // find the index of this pixel in the mask vector
                int mask_index = posToIndex(x_iter, y_iter, mask_diameter);
                // find the index of this pixel in the canvas vector
                int canvas_index = posToIndex(canvas_x, canvas_y, m_width);

                // retrieve canvas colour
                RGBA canvas_colour = m_data[canvas_index];

                // retrieve mask and canvas intensities, and brush colour
                auto [mask_intensity, canvas_intensity, curr_colour] = getIntensitiesAndColour(mask_index, canvas_index);

                uint8_t new_red = mask_intensity * curr_colour.r + canvas_intensity * canvas_colour.r;
                uint8_t new_green = mask_intensity * curr_colour.g + canvas_intensity * canvas_colour.g;
                uint8_t new_blue = mask_intensity * curr_colour.b + canvas_intensity * canvas_colour.b;

                m_data[canvas_index] = RGBA{new_red, new_green, new_blue};
            }

        }
    }
}

/**
 * @brief The following functions give the mask intensities for a pixel of a given
 *        distance from the center of the mask.
 */
float constIntensity(float dist) {
    if (dist > settings.brushRadius) {
        return 0;
    } else {
        return 1;
    }
}

float linIntensity(float dist) {
    if (dist > settings.brushRadius) {
        return 0;
    } else {
        return 1 - dist / (float)settings.brushRadius;
    }
}

float quadIntensity(float dist) {
    int r = settings.brushRadius;
    if (dist > r) {
        return 0;
    } else {
        float a = (float)1 / (float)pow(r, 2);
        float b= (float)2 / (float)r;
        return (a*pow(dist, 2)) - (b*dist) + 1;
    }
}

/**
 * @brief Returns a 1D vector which represents the intensity values in the mask.
 *        Generated according to current brush settings.
 */
void Canvas2D::changeMask(float (&intensity) (float)) {
    int mask_diameter = 2 * settings.brushRadius + 1;

    // clear current mask
    mask_vector.clear();
    // reserve (2R+1)^2 elements in the mask vector
    mask_vector.reserve(pow(mask_diameter, 2));

    // populate the mask vector
    for (int y = 0; y < mask_diameter; y++) {
        for (int x = 0; x < mask_diameter; x++) {
            mask_vector.push_back(intensity(maskDist(x, y)));
        }
    }
}

/**
 * @brief based on the current settings.brushType value, changes the mask accordingly.
 */
void Canvas2D::setMask() {
    switch(settings.brushType) {
    case 0: // constant
        changeMask(constIntensity);
        break;
    case 1: // linear
        changeMask(linIntensity);
        break;
    case 2: // quadratic
        changeMask(quadIntensity);
        break;
    case 3: // smudge
        // arbitrarily chose linear mask for smudge, can be anything really
        changeMask(linIntensity);
        break;
    case 4: // spray
        // choose constant intensity
        changeMask(constIntensity);
        break;
    default:
        break;
    }
}

/**
 * @brief stores the pixels around the current position in the smudge_data
 *        vector, to be placed down when the mouse is dragged.
 */
void Canvas2D::storeSmudgeData(int x, int y) {
    int r = settings.brushRadius;
    int mask_diameter = 2 * r + 1;

    // clear current smudge data
    smudge_data.clear();
    // reserve (2R+1)^2 elements in the smudge data vector
    smudge_data.reserve(pow(mask_diameter, 2));

    for (int y_iter = y - r; y_iter < y - r + mask_diameter; y_iter++) {
        for (int x_iter = x - r; x_iter < x - r + mask_diameter; x_iter++) {
            // store the current canvas pixel
            smudge_data.push_back(m_data[posToIndex(x_iter, y_iter, m_width)]);
        }
    }
}

/**
 * @brief These functions are called when the mouse is clicked and dragged on the canvas
 */
void Canvas2D::mouseDown(int x, int y) {
    if (settings.brushType == 3) { // if smudge brush
        // store smudge data
        storeSmudgeData(x, y);
    }
    paintBrush(x, y);
    displayImage();
}

void Canvas2D::mouseDragged(int x, int y) {
    paintBrush(x, y);
    if (settings.brushType == 3) { // if smudge brush
        // store smudge data
        storeSmudgeData(x, y);
    }
    displayImage();
}

void Canvas2D::mouseUp(int x, int y) {

}
