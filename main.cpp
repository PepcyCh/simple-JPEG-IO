#include "JPEGUtility.hpp"
#include <iostream>

int main() {
    Image img = JPEGUtility::JPEG()->read("rgb.jpg");

    std::cout << "width = " << img.getWidth() << ", height = " << img.getHeight() << std::endl;

    JPEGUtility::JPEG()->write(img, "copy.jpg");
    
    return 0;
}
