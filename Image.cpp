#include "Image.hpp"

#include <algorithm>

Image::Image(size_t width, size_t height) : width(width), height(height) {
    data.resize(width * height * 4);
}

Image::~Image() {
    data.clear();
}

Image::Image(const Image& rhs) {
    width = rhs.width;
    height = rhs.height;
    data = rhs.data;
}

Image& Image::operator=(const Image& rhs) {
    width = rhs.width;
    height = rhs.height;
    data = rhs.data;
    return *this;
}

bool Image::empty() const { return data.empty(); }

void Image::resize(size_t width, size_t height) {
    data.clear();
    data.resize(width * height * 4);
    this->width = width;
    this->height = height;
}

size_t Image::getWidth() const { return width; }
size_t Image::getHeight() const { return height; }

uint8_t Image::get(int x, int y, int i) const {
    x = std::clamp<int>(x, 0, height - 1);
    y = std::clamp<int>(y, 0, width - 1);
    return data[(x * width + y) * 4 + i];
}

uint8_t Image::getR(int x, int y) const {
    x = std::clamp<int>(x, 0, height - 1);
    y = std::clamp<int>(y, 0, width - 1);
    return data[(x * width + y) * 4];
}
uint8_t Image::getG(int x, int y) const {
    x = std::clamp<int>(x, 0, height - 1);
    y = std::clamp<int>(y, 0, width - 1);
    return data[(x * width + y) * 4 + 1];
}
uint8_t Image::getB(int x, int y) const {
    x = std::clamp<int>(x, 0, height - 1);
    y = std::clamp<int>(y, 0, width - 1);
    return data[(x * width + y) * 4 + 2];
}
uint8_t Image::getA(int x, int y) const {
    x = std::clamp<int>(x, 0, height - 1);
    y = std::clamp<int>(y, 0, width - 1);
    return data[(x * width + y) * 4 + 3];
}

void Image::setColor(int x, int y, int i, uint8_t v) {
    if (x < 0 || y < 0 || x >= height || y >= width) return;
    data[(x * width + y) * 4 + i] = v;
}

void Image::setColor(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
    if (x < 0 || y < 0 || x >= height || y >= width) return;
    size_t id = (x * width + y) * 4;
    data[id] = r;
    data[id + 1] = g;
    data[id + 2] = b;
    data[id + 3] = 255;
}

void Image::setColor(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
    if (x < 0 || y < 0 || x >= height || y >= width) return;
    size_t id = (x * width + y) * 4;
    data[id] = r;
    data[id + 1] = g;
    data[id + 2] = b;
    data[id + 3] = a;
}
