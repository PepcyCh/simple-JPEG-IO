#pragma once

#include <cstdlib>
#include <vector>

class Image {
public:
    Image(size_t width, size_t height);
    ~Image();

    Image(const Image& rhs);
    Image& operator=(const Image& rhs);

    bool empty() const;

    void resize(size_t width, size_t height);

    size_t getWidth() const;
    size_t getHeight() const;

    uint8_t get(int x, int y, int i) const;

    uint8_t getR(int x, int y) const;
    uint8_t getG(int x, int y) const;
    uint8_t getB(int x, int y) const;
    uint8_t getA(int x, int y) const;

    void setColor(int x, int y, int i, uint8_t v);
    void setColor(int x, int y, uint8_t r, uint8_t g, uint8_t b);
    void setColor(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a);

private:
    size_t width;
    size_t height;
    std::vector<uint8_t> data;
};
