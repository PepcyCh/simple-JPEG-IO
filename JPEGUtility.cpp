#include "JPEGUtility.hpp"

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <queue>
#include "Mat.hpp"

// #define JPEG_DUMP

size_t HuffmanTable::buildHuffmanTableFromDHT(const uint8_t* num, const uint8_t* data) {
    std::queue<std::shared_ptr<HuffmanNode>> nodes;

    size_t total = 0;
    for (int i = 0; i < 16; i++) total += num[i];

    size_t offset = total;
    for (int l = 15; l >= 0; l--) {
        offset -= num[l];
        size_t childrenCount = nodes.size();

        for (int i = 0; i < num[l]; i++) {
            nodes.push(std::make_shared<HuffmanNode>(data[offset + i]));
        }

        while (childrenCount > 0) {
            std::shared_ptr<HuffmanNode> node = std::make_shared<HuffmanNode>();
            std::shared_ptr<HuffmanNode> lc = nodes.front();
            nodes.pop();
            --childrenCount;
            node->c[0] = lc;
            if (childrenCount > 0) {
                std::shared_ptr<HuffmanNode> rc = nodes.front();
                nodes.pop();
                --childrenCount;
                node->c[1] = rc;
            }
            nodes.push(node);
        }
    }

    root = std::make_shared<HuffmanNode>();
    std::shared_ptr<HuffmanNode> lc = nodes.front();
    nodes.pop();
    root->c[0] = lc;
    if (!nodes.empty()) {
        std::shared_ptr<HuffmanNode> rc = nodes.front();
        nodes.pop();
        root->c[1] = rc;
    }

    return total;
}

void HuffmanTable::buildSimpleDCHuffmanTable(uint8_t* num, std::vector<uint8_t>& data, std::pair<int, int>* codeMap) {
    std::fill(num, num + 16, 0);
    data.clear();
    for (int i = 0; i < 12; i++) {
        data.push_back(i);
        ++num[3];
    }
    buildHuffmanTableFromDHT(num, data.data());
    getCodeMap(codeMap);
}
void HuffmanTable::buildSimpleACHuffmanTable(uint8_t* num, std::vector<uint8_t>& data, std::pair<int, int>* codeMap) {
    std::fill(num, num + 16, 0);
    data.clear();
    data.push_back(0);
    ++num[7];
    data.push_back(0xF0);
    ++num[7];
    for (int r = 0; r < 16; r++) {
        for (int s = 1; s <= 10; s++) {
            int rs = (r << 4) | s;
            data.push_back(rs);
            ++num[7];
        }
    }
    buildHuffmanTableFromDHT(num, data.data());
    getCodeMap(codeMap);
}

uint16_t HuffmanTable::parse(const uint8_t* data, size_t size, size_t& byteOffset, size_t& bitOffset) const {
    std::shared_ptr<HuffmanNode> current = root;
    for (size_t i = byteOffset; i < size; i++) {
        for (int j = 7 - bitOffset; j >= 0; j--) {
            size_t x = (data[i] >> j) & 1u;
            current = current->c[x];
            assert(current);
            if (!current) {
                byteOffset = size;
                return -1;
            }
            if (current->isLeaf) {
                if (j == 0) {
                    byteOffset = i + 1;
                    bitOffset = 0;
                } else {
                    byteOffset = i;
                    bitOffset = 8 - j;
                }
                return current->value;
            }
        }
        bitOffset = 0;
    }
    return 0;
}

void HuffmanTable::getCodeMap(std::pair<int, int>* codeMap) {
    std::fill(codeMap, codeMap + 256, std::make_pair(0, 0));
    std::queue<std::pair<std::shared_ptr<HuffmanNode>, std::pair<int, int>>> q;
    q.emplace(root, std::make_pair(0, 0));
    while (!q.empty()) {
        auto u = q.front();
        q.pop();
        if (u.first->isLeaf) {
            codeMap[u.first->value] = u.second;
        } else {
            if (u.first->c[0])
                q.emplace(u.first->c[0], std::make_pair(u.second.first + 1, u.second.second << 1));
            if (u.first->c[1])
                q.emplace(u.first->c[1], std::make_pair(u.second.first + 1, u.second.second << 1 | 1));
        }
    }
}

// ------------------

static void YCbCr2RGB(uint8_t Y, uint8_t Cb, uint8_t Cr, uint8_t& R, uint8_t& G, uint8_t& B) {
    R = std::clamp(Y + 1.402 * (Cr - 128), 0., 255.);
    G = std::clamp(Y - 0.344136 * (Cb - 128) - 0.714136 * (Cr - 128), 0., 255.);
    B = std::clamp(Y + 1.772 * (Cb - 128), 0., 255.);
}
static void RGB2YCbCr(uint8_t R, uint8_t G, uint8_t B, uint8_t& Y, uint8_t& Cb, uint8_t& Cr) {
    Y = 0.299 * R + 0.587 * G + 0.114 * B;
    Cb = 128 - 0.168736 * R - 0.331264 * G + 0.5 * B;
    Cr = 128 + 0.5 * R - 0.418688 * G - 0.081312 * B;
}

static uint16_t parseBigendian16(uint16_t x) { return ((x >> 8u) | ((x & 0xFFu) << 8u)); }
static uint16_t getBigendian16(uint16_t x) { return ((x >> 8u) | ((x & 0xFFu) << 8u)); }
static const char* getBigendian16PChar(uint16_t x) {
    static char c[2];
    c[0] = (x >> 8) & 0xFF;
    c[1] = x & 0xFF;
    return c;
}

static Mat<int> DCT(const Mat<int>& mat) {
    /*
    Mat<int> res(8, 8);
    for (int u = 0; u < 8; u++) {
        double cu = u ? 1 : 1.0 / sqrt(2);
        for (int v = 0; v < 8; v++) {
            double cv = v ? 1 : 1.0 / sqrt(2);
            double sum = 0;
            for (int x = 0; x < 8; x++) {
                for (int y = 0; y < 8; y++) {
                    sum += mat[y][x] *
                           cos((2 * x + 1) * u * M_PI / 16.0) *
                           cos((2 * y + 1) * v * M_PI / 16.0);
                }
            }
            sum = cu * cv * sum / 4.0;
            res[v][u] = round(sum);
        }
    }
    return res;
     */
    Mat<double> temp(8, 8, 0);
    for (int u = 0; u < 8; u++) {
        for (int v = 0; v < 8; v++) {
            double sum = 0;
            for (int i = 0; i < 8; i++) {
                sum += mat[i][v] * cos((2 * i + 1) * u * M_PI / 16.0);
            }
            temp[u][v] = sum;
        }
    }
    Mat<int> res(8, 8, 0);
    for (int v = 0; v < 8; v++) {
        double cv = v ? 1 : 1.0 / sqrt(2);
        for (int u = 0; u < 8; u++) {
            double cu = u ? 1 : 1.0 / sqrt(2);
            double sum = 0;
            for (int i = 0; i < 8; i++) {
                sum += temp[u][i] * cos((2 * i + 1) * v * M_PI / 16.0);
            }
            sum = cu * cv * sum / 4.0;
            res[u][v] = round(sum);
        }
    }
    return res;
}

static Mat<int> IDCT(const Mat<int>& mat) {
    /*
    Mat<int> res(8, 8, 0);
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            double sum = 0;
            for (int u = 0; u < 8; u++) {
                double cu = u ? 1 : 1.0 / sqrt(2);
                for (int v = 0; v < 8; v++) {
                    double cv = v ? 1 : 1.0 / sqrt(2);
                    sum += cu * cv * mat[v][u] *
                           cos((2 * x + 1) * u * M_PI / 16.0) *
                           cos((2 * y + 1) * v * M_PI / 16.0);
                }
            }
            sum /= 4.0;
            res[y][x] = round(sum);
        }
    }
    return res;
     */
    Mat<double> temp(8, 8, 0);
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            double sum = 0;
            for (int i = 0; i < 8; i++) {
                double ci = i ? 1 : 1.0 / sqrt(2);
                sum += ci * mat[x][i] * cos((2 * y + 1) * i * M_PI / 16.0);
            }
            temp[x][y] = sum;
        }
    }
    Mat<int> res(8, 8, 0);
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            double sum = 0;
            for (int i = 0; i < 8; i++) {
                double ci = i ? 1 : 1.0 / sqrt(2);
                sum += ci * temp[i][y] * cos((2 * x + 1) * i * M_PI / 16.0);
            }
            sum /= 4.0;
            res[x][y] = round(sum);
        }
    }
    return res;
}

// --------------------

class Component {
public:
    ~Component() {
        delete[] data;
        delete[] blocks;
    }

    void resize(int width, int height, bool isProgress = false) {
        this->width = width;
        this->height = height;
        delete[] data;
        delete[] blocks;
        data = new uint8_t[width * height];
        blocks = nullptr;
        _mcuX = (width + 7) >> 3;
        _mcuY = (height + 7) >> 3;
        _mcuCount = _mcuX * _mcuY;
        mcuIndex = mcuIndexTemp = 0;
        if (isProgress) {
            blocks = new Mat<int>[_mcuCount];
            for (int i = 0; i < _mcuCount; i++) blocks[i].create(8, 8, 0);
        }
    }

    void setAt(int x, int y, int value) {
        if (x < 0 || x >= height || y < 0 || y >= width) return;
        data[x * width + y] = std::clamp(value, 0, 255);
    }

    void setBlockAt(int mcux, int mcuy, int i, int j, int value) {
        if (mcux < 0 || mcux >= _mcuX || mcuy < 0 || mcuy >= _mcuY) return;
        blocks[mcuy * _mcuX + mcux][i][j] = value;
    }
    void setBlockAt(int mcuid, int i, int j, int value) {
        blocks[mcuid][i][j] = value;
    }
    void addBlockAt(int mcux, int mcuy, int i, int j, int value) {
        if (mcux < 0 || mcux >= _mcuX || mcuy < 0 || mcuy >= _mcuY) return;
        blocks[mcuy * _mcuX + mcux][i][j] += value;
    }
    void addBlockAt(int mcuid, int i, int j, int value) {
        blocks[mcuid][i][j] += value;
    }
    int getBlockAt(int mcux, int mcuy, int i, int j) const {
        if (mcux < 0 || mcux >= _mcuX || mcuy < 0 || mcuy >= _mcuY) return 0;
        return blocks[mcuy * _mcuX + mcux][i][j];
    }
    int getBlockAt(int mcuid, int i, int j) const {
        return blocks[mcuid][i][j];
    }

    void DequantAndIDCT(const Mat<uint16_t>& quant) {
        for (int i = 0; i < _mcuCount; i++) {
            DequantndIDCTBlock(i, quant);
        }
    }
    void DCTAndQuant(const Mat<uint16_t>& quant) {
        for (int i = 0; i < _mcuCount; i++) {
            DCTAndQuantBlock(i, quant);
        }
    }

    uint8_t sampleAt(double u, double v) {
        float tu = std::clamp(u, 0.0, 0.99999) * height;
        float tv = std::clamp(v, 0.0, 0.99999) * width;

        /*
        // bilinear
        int su[2];
        su[0] = std::clamp<int>(round(tu) - 1, 0, height - 1);
        su[1] = std::clamp<int>(su[0] + 1, 0, height - 1);
        float du = tu - 0.5f - su[0];
        if (du < 0) su[1] = su[0];

        int sv[2];
        sv[0] = std::clamp<int>(round(tv) - 1, 0, width - 1);
        sv[1] = std::clamp<int>(sv[0] + 1, 0, width - 1);
        float dv = tv - 0.5f - sv[0];
        if (dv < 0) sv[1] = sv[0];

        double sum = 0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                uint8_t t = data[su[i] * width + sv[j]];
                sum += t;
            }
        }
        return (uint8_t) std::clamp(round(sum / 4), 0., 255.);
         */
        // nearest
        int x = floor(tu);
        int y = floor(tv);
        return (uint8_t) data[x * width + y];
    }

    int hFactor() const { return horizontalSamplingFactor; }
    void setHFactor(int w) { horizontalSamplingFactor = w; }
    int vFactor() const { return verticalSamplingFactor; }
    void setVFactor(int w) { verticalSamplingFactor = w; }
    int quantization() const { return quantizationSelector; }
    void setQuantizationSelector(int quantization) { quantizationSelector = quantization; }
    int DCHuffman() const { return DCHuffmanSelector; }
    void setDChuffmanSelector(int DCHuffman) { DCHuffmanSelector = DCHuffman; }
    int ACHuffman() const { return ACHuffmanSelector; }
    void setAChuffmanSelector(int ACHuffman) { ACHuffmanSelector = 4 | ACHuffman; }
    int getWidth() const { return width; }
    int getHeight() const { return height; }

    int mcuIndex;
    int mcuIndexTemp;

    int mcuX() const { return _mcuX; }
    int mcuY() const { return _mcuY; }
    int mcuCount() const { return _mcuCount; }

private:
    int horizontalSamplingFactor;
    int verticalSamplingFactor;
    int quantizationSelector;
    int DCHuffmanSelector;
    int ACHuffmanSelector;

    int width, height;
    int _mcuX, _mcuY, _mcuCount;
    uint8_t* data = nullptr;
    Mat<int>* blocks = nullptr;

    void DequantndIDCTBlock(int mcuid, const Mat<uint16_t>& quant) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                blocks[mcuid][i][j] *= quant[i][j];
            }
        }
        Mat<int> temp = IDCT(blocks[mcuid]);
        int x = mcuid / _mcuX;
        int y = mcuid % _mcuX;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                setAt(x * 8 + i, y * 8 + j, temp[i][j] + 128);
            }
        }
    }
    void DCTAndQuantBlock(int mcuid, const Mat<uint16_t>& quant) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                blocks[mcuid][i][j] -= 128;
            }
        }
        Mat<int> temp = DCT(blocks[mcuid]);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                blocks[mcuid][i][j] = temp[i][j] / quant[i][j];
            }
        }
    }
};
static Component compos[4];

// --------------------

JPEGUtility* JPEGUtility::singleton = nullptr;

void JPEGUtility::reconstruct(Image& img) {
    if (isProgress) {
        for (int cid = 0; cid < componentsCount; cid++) {
            compos[cid].DequantAndIDCT(quantizations[compos[cid].quantization()]);
        }
    }

    for (int i = 0; i < img.getHeight(); i++) {
        for (int j = 0; j < img.getWidth(); j++) {
            double u = (double) i / img.getHeight();
            double v = (double) j / img.getWidth();
            uint8_t Y = compos[0].sampleAt(u, v);
            uint8_t Cb = compos[1].sampleAt(u, v);
            uint8_t Cr = compos[2].sampleAt(u, v);
            uint8_t A = componentsCount > 3 ? compos[3].sampleAt(u, v) : 255;
            uint8_t R, G, B;
            YCbCr2RGB(Y, Cb, Cr, R, G, B);
            img.setColor(i, j, R, G, B, A);
        }
    }
}

int JPEGUtility::getExtendValue(size_t& byteOffset, size_t& bitOffset, size_t length) {
    if (length == 0) return 0;
    int temp;
    if (bitOffset + length <= 8) {
        temp = (scanData[byteOffset] >> (8 - bitOffset - length) & ((1 << length) - 1));
    } else {
        size_t bitsInFirstByte = 8 - bitOffset;
        size_t fullBytes = (length - bitsInFirstByte) / 8;
        size_t bitsInLastByte = length - bitsInFirstByte - 8 * fullBytes;

        temp = (scanData[byteOffset] & ((1 << bitsInFirstByte) - 1));
        for (int i = 1; i <= fullBytes; i++) {
            temp = (temp << 8) | scanData[byteOffset + i];
        }
        temp = (temp << bitsInLastByte) |
               (scanData[byteOffset + fullBytes + 1] >> (8 - bitsInLastByte));
    }

    int value;
    if (length > 0 && ((temp >> (length - 1)) & 1) == 0) {
        value = -(int) (~temp & ((1 << length) - 1));
    } else {
        value = temp;
    }

    bitOffset += length;
    while (bitOffset >= 8) {
        ++byteOffset;
        bitOffset -= 8;
    }

    return value;
}

void JPEGUtility::pushValue(size_t& byteOffset, size_t& bitOffset, const std::pair<int, int>& code) {
    if (code.first == 0) return;
    int newByteOffset = byteOffset, newBitOffset = bitOffset + code.first;
    while (newBitOffset >= 8) {
        ++newByteOffset;
        newBitOffset -= 8;
    }
    while (scanData.size() <= newByteOffset) scanData.push_back(0);
    if (bitOffset + code.first <= 8) {
        scanData[byteOffset] |= (code.second << (8 - bitOffset - code.first));
    } else {
        size_t bitsInFirstByte = 8 - bitOffset;
        size_t fullBytes = (code.first - bitsInFirstByte) / 8;
        size_t bitsInLastByte = code.first - bitsInFirstByte - 8 * fullBytes;
        scanData[byteOffset] |= (code.second >> (code.first - bitsInFirstByte));
        for (int i = 1; i <= fullBytes; i++) {
            scanData[byteOffset + i] = (code.second >> (code.first - bitsInFirstByte - 8 * i)) & 0xFF;
        }
        scanData[byteOffset + fullBytes + 1] = (code.second << (8 - bitsInLastByte)) & 0xFF;
    }
    byteOffset = newByteOffset;
    bitOffset = newBitOffset;
}

int JPEGUtility::getCodeLength(int val) {
    if (val < 0) val = -val;
    assert(val < (1 << 15));
    for (int i = 15; i >= 0; i--) if (val & (1 << i)) return i + 1;
    return 0;
}

int JPEGUtility::getExtendCode(int val) {
    if (val >= 0) {
        return val;
    } else {
        val = -val;
        int l = getCodeLength(val);
        return (~val) & ((1 << l) - 1);
    }
}

Mat<int> JPEGUtility::parseECSBlock_base(int cid, size_t& byteOffset, size_t& bitOffset) {
    Mat<int> block(8, 8, 0);

    // parse DC
    {
        uint16_t value = huffmans[compos[cid].DCHuffman()].parse(
                scanData.data(), scanData.size(), byteOffset, bitOffset);
        uint16_t DCLength = value & 0x0F;

        int DC = getExtendValue(byteOffset, bitOffset, DCLength);

        DC += predDC[cid];
        predDC[cid] = DC;

        block[0][0] = DC;
    }

    // parse AC
    size_t ACIndex = 1;
    while (byteOffset < scanData.size() && ACIndex < 64) {
        uint16_t value = huffmans[compos[cid].ACHuffman()].parse(
                scanData.data(), scanData.size(), byteOffset, bitOffset);

        if (value == 0) { // EOB
            break;
        } else if (value == 0xF0) { // ZRL
            ACIndex += 16;
            continue;
        }
        ACIndex += (value >> 4);

        uint16_t ACLength = value & 0x0F;
        int AC = getExtendValue(byteOffset, bitOffset, ACLength);

        size_t aid = zz2quad[ACIndex];
        block[aid >> 3][aid & 0x07] = AC;
        ++ACIndex;
    }

    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            block[x][y] *= quantizations[compos[cid].quantization()][x][y];
        }
    }
    Mat<int> ret = IDCT(block);
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            block[x][y] = ret[x][y] + 128;
        }
    }

    return block;
}

void JPEGUtility::parseECS_base() {
    std::fill(predDC, predDC + 4, 0);
    size_t byteOffset = 0, bitOffset = 0;
    while (byteOffset < scanData.size() && mcuIndex < mcuCount) {
        for (int i = 0; i < scanCompoTable.size(); i++) {
            int cid = compoMap[scanCompoTable[i].scanComponentsSelector];

            for (int u = 0; u < compos[cid].vFactor(); u++) {
                for (int v = 0; v < compos[cid].hFactor(); v++) {
                    Mat<int> block = parseECSBlock_base(cid, byteOffset, bitOffset);

                    size_t mcuX = mcuIndex % mcuCountX;
                    size_t mcuY = mcuIndex / mcuCountX;
                    size_t subX = (mcuX * compos[cid].hFactor() + v) * 8;
                    size_t subY = (mcuY * compos[cid].vFactor() + u) * 8;
                    for (int x = 0; x < 8; x++) {
                        for (int y = 0; y < 8; y++) {
                            compos[cid].setAt(subY + x, subX + y, block[x][y]);
                        }
                    }
                }
            }
        }

        ++mcuIndex;
        if (restartInterval != 0 && mcuIndex % restartInterval == 0) {
            if (bitOffset) {
                bitOffset = 0;
                ++byteOffset;
            }
            assert(byteOffset == scanData.size());
            break;
        }
    }
}

void JPEGUtility::parseECS_prog_DC() {
    size_t byteOffset = 0, bitOffset = 0;
    if (pSOSEnd->successiveApproximationBitPositionHigh() == 0) { // first
        std::fill(predDC, predDC + 4, 0);
        while (byteOffset < scanData.size() && mcuIndex < mcuCount) {
            for (int i = 0; i < scanCompoTable.size(); i++) {
                int cid = compoMap[scanCompoTable[i].scanComponentsSelector];

                for (int u = 0; u < compos[cid].vFactor(); u++) {
                    for (int v = 0; v < compos[cid].hFactor(); v++) {
                        uint16_t value = huffmans[compos[cid].DCHuffman()].parse(
                                scanData.data(), scanData.size(), byteOffset, bitOffset);
                        uint16_t DCLength = value & 0x0F;

                        int DC = getExtendValue(byteOffset, bitOffset, DCLength);

                        DC += predDC[cid];
                        predDC[cid] = DC;

                        size_t mcuX = mcuIndex % mcuCountX;
                        size_t mcuY = mcuIndex / mcuCountX;
                        size_t subX = mcuX * compos[cid].hFactor() + v;
                        size_t subY = mcuY * compos[cid].vFactor() + u;
                        compos[cid].setBlockAt(subX, subY, 0, 0,
                                DC << pSOSEnd->successiveApproximationBitPositionLow());
                    }
                }
            }

            ++mcuIndex;
            if (restartInterval != 0 && mcuIndex % restartInterval == 0) {
                if (bitOffset) {
                    bitOffset = 0;
                    ++byteOffset;
                }
                assert(byteOffset == scanData.size());
                break;
            }
        }
    } else {
        int refine = 1 << pSOSEnd->successiveApproximationBitPositionLow();
        while (byteOffset < scanData.size() && mcuIndex < mcuCount) {
            for (int i = 0; i < scanCompoTable.size(); i++) {
                int cid = compoMap[scanCompoTable[i].scanComponentsSelector];

                for (int u = 0; u < compos[cid].vFactor(); u++) {
                    for (int v = 0; v < compos[cid].hFactor(); v++) {
                        size_t mcuX = mcuIndex % mcuCountX;
                        size_t mcuY = mcuIndex / mcuCountX;
                        size_t subX = mcuX * compos[cid].hFactor() + v;
                        size_t subY = mcuY * compos[cid].vFactor() + u;

                        int value = getExtendValue(byteOffset, bitOffset);
                        if (value > 0) compos[cid].addBlockAt(subX, subY, 0, 0,
                                               compos[cid].getBlockAt(subX, subY, 0, 0) > 0 ? refine : -refine);
                    }
                }
            }

            ++mcuIndex;
            if (restartInterval != 0 && mcuIndex % restartInterval == 0) {
                if (bitOffset) {
                    bitOffset = 0;
                    ++byteOffset;
                }
                assert(byteOffset == scanData.size());
                break;
            }
        }
    }
}

void JPEGUtility::parseECS_prog_AC() {
    size_t byteOffset = 0, bitOffset = 0;
    int cid = compoMap[scanCompoTable[0].scanComponentsSelector];
    compos[cid].mcuIndex = compos[cid].mcuIndexTemp;
    int start = pSOSEnd->startOfSpectralOfPredictorSelection;
    int end = pSOSEnd->endOfSpectralOfPredictorSelection;

    if (pSOSEnd->successiveApproximationBitPositionHigh() == 0) { // first
        int EOBs = 0;
        while (byteOffset < scanData.size() && compos[cid].mcuIndex < compos[cid].mcuCount()) {
            size_t ACIndex = start;

            while (byteOffset < scanData.size() && ACIndex <= end) {
                if (EOBs > 0) {
                    --EOBs;
                    break;
                }

                uint16_t value = huffmans[compos[cid].ACHuffman()].parse(
                        scanData.data(), scanData.size(), byteOffset, bitOffset);
                uint16_t ACLength = value & 0x0F;

                if (!ACLength) {
                    uint16_t r = value >> 4;
                    if (r < 15) { // EOB
                        EOBs = getExtendValue(byteOffset, bitOffset, r);
                        if (EOBs >= 0) EOBs += (1 << r);
                        else EOBs += ((1 << (r + 1)) - 1);
                        --EOBs;
                        break;
                    }
                    // ZRL
                    ACIndex += 16;
                    continue;
                }

                ACIndex += (value >> 4);

                int AC = getExtendValue(byteOffset, bitOffset, ACLength);

                size_t aid = zz2quad[ACIndex];
                compos[cid].setBlockAt(compos[cid].mcuIndex, aid >> 3, aid & 0x07,
                        AC << pSOSEnd->successiveApproximationBitPositionLow());
                ++ACIndex;
            }

            ++compos[cid].mcuIndex;
            if (restartInterval != 0 && compos[cid].mcuIndex % restartInterval == 0) {
                if (bitOffset) {
                    bitOffset = 0;
                    ++byteOffset;
                }
                assert(byteOffset == scanData.size());
                break;
            }
        }
    } else {
        int EOBs = 0;
        int refine = 1 << pSOSEnd->successiveApproximationBitPositionLow();
        while (byteOffset < scanData.size() && compos[cid].mcuIndex < compos[cid].mcuCount()) {
            size_t ACIndex = start;

            while (byteOffset < scanData.size() && ACIndex <= end) {
                if (EOBs > 0) {
                    --EOBs;
                    for (size_t k = start; k <= end; k++) {
                        int his, aid = zz2quad[k];
                        if ((his = compos[cid].getBlockAt(compos[cid].mcuIndex, aid >> 3, aid & 0x07))) {
                            int v = getExtendValue(byteOffset, bitOffset);
                            if (v > 0 && !(his & refine)) compos[cid].addBlockAt(compos[cid].mcuIndex, aid >> 3, aid & 0x07,
                                    his > 0 ? refine : -refine);
                        }
                    }
                    break;
                }

                uint16_t value = huffmans[compos[cid].ACHuffman()].parse(
                        scanData.data(), scanData.size(), byteOffset, bitOffset);
                uint16_t length = value & 0x0F;

                int r = value >> 4, AC = 0;
                if (!length) {
                    if (r < 15) { // EOBs
                        EOBs = getExtendValue(byteOffset, bitOffset, r);
                        if (EOBs >= 0) EOBs += (1 << r);
                        else EOBs += ((1 << (r + 1)) - 1);
                        --EOBs;
                        r = 64;
                    }
                } else {
                    int v = getExtendValue(byteOffset, bitOffset);
                    AC = v > 0 ? refine : -refine;
                }
                for (; byteOffset < scanData.size() && ACIndex <= end; ACIndex++) {
                    int his, aid = zz2quad[ACIndex];
                    if ((his = compos[cid].getBlockAt(compos[cid].mcuIndex, aid >> 3, aid & 0x07))) {
                        int v = getExtendValue(byteOffset, bitOffset);
                        if (v > 0 && !(his & refine)) compos[cid].addBlockAt(compos[cid].mcuIndex, aid >> 3, aid & 0x07,
                                his > 0 ? refine : -refine);
                    } else {
                        if (r == 0) {
                            compos[cid].setBlockAt(compos[cid].mcuIndex, aid >> 3, aid & 0x07, AC);
                            ++ACIndex;
                            break;
                        } else {
                            --r;
                        }
                    }
                }
            }

            ++compos[cid].mcuIndex;
            if (restartInterval != 0 && compos[cid].mcuIndex % restartInterval == 0) {
                if (bitOffset) {
                    bitOffset = 0;
                    ++byteOffset;
                }
                assert(byteOffset == scanData.size());
                break;
            }
        }
    }
}

size_t JPEGUtility::parseECS(uint8_t* data, uint8_t* end) {
    scanData.clear();
    size_t scanLength = 0;

    {
        uint8_t *p = data;
        bool bitstuff = false;
        while (p < end && (*p != 0xFF || *(p + 1) == 0x00)) {
            if (bitstuff) {
                bitstuff = false;
            } else {
                scanData.push_back(*p);
            }
            bitstuff = (*p == 0xFF && *(p + 1) == 0x00);

            ++p;
            ++scanLength;
        }
    }

    mcuIndex = mcuIndexTemp;
    if (!isProgress) {
        parseECS_base();
    } else if (pSOSEnd->startOfSpectralOfPredictorSelection == 0) {
        parseECS_prog_DC();
    } else {
        parseECS_prog_AC();
    }

    return scanLength;
}

Image JPEGUtility::read(const std::string &filename) {
    std::ifstream fp(filename, std::ios::binary);
    if (!fp.is_open()) {
        std::cerr << "Failed to open file '" << filename << "'" << std::endl;
        return Image(0, 0);
    }

    fp.seekg(0, std::ios::end);
    size_t filesize = fp.tellg();
    fp.seekg(0);

    uint8_t* filedata = new uint8_t[filesize];
    fp.read(reinterpret_cast<char*>(filedata), filesize);
    fp.close();

    JPEGSOI* pSOI = reinterpret_cast<JPEGSOI*>(filedata);
    if (pSOI->SOI != 0xD8FF) {
        std::cerr << "Not a JPEG file" << std::endl;
        delete[] filedata;
        return Image(0, 0);
    }

    uint8_t* pData = filedata + 2;
    uint8_t* pEnd = filedata + filesize;
    restartInterval = 0;

    Image img(0, 0);
    while (pData < pEnd) {
        JPEGMarkerHeader* pHeader = reinterpret_cast<JPEGMarkerHeader*>(pData);
        switch (pHeader->marker) {
            case 0xD9FF: {
#ifdef JPEG_DUMP
                std::cerr << "End of Image" << std::endl;
#endif
                goto endOfParse;
            }
            case 0xC0FF:
            case 0xC2FF: {
                if (pHeader->marker == 0xC0FF) {
#ifdef JPEG_DUMP
                    std::cerr << "Start of Frame (baseline DCT)" << std::endl;
                    std::cerr << "-----------------------------" << std::endl;
#endif
                    isProgress = false;
                } else {
#ifdef JPEG_DUMP
                    std::cerr << "Start of Frame (progressive DCT)" << std::endl;
                    std::cerr << "--------------------------------" << std::endl;
#endif
                    isProgress = true;
                }

                JPEGSOF* pSOF = reinterpret_cast<JPEGSOF*>(pData);

                samplePrecision = pSOF->samplePrecision;
                width = parseBigendian16(pSOF->numberOfSamplesPerLine);
                height = parseBigendian16(pSOF->numberOfLines);
                img.resize(width, height);
                mcuCountX = (width + 7) >> 3;
                mcuCountY = (height + 7) >> 3;
                mcuIndex = mcuIndexTemp = 0;
                componentsCount = pSOF->numberOfComponents;
#ifdef JPEG_DUMP
                std::cerr << "Sample Precision: " << (int) pSOF->samplePrecision << std::endl;
                std::cerr << "Number of Lines: " << (int) height << std::endl;
                std::cerr << "Number of Samples per Line: " << (int) width << std::endl;
                std::cerr << "Number of Components: " << (int) pSOF->numberOfComponents << std::endl;
#endif
                frameCompoTable.clear();
                uint8_t* temp = pData + sizeof(JPEGSOF);
                int maxHFactor = 1, maxVFactor = 1;
                for (int i = 0; i < pSOF->numberOfComponents; i++) {
                    JPEGFrameSpecification* pCompo = reinterpret_cast<JPEGFrameSpecification*>(temp);
                    frameCompoTable.push_back(*pCompo);
                    temp += sizeof(JPEGFrameSpecification);
#ifdef JPEG_DUMP
                    std::cerr << "  Component Identifier: " << (int) pCompo->componentsIdentifier << std::endl;
                    std::cerr << "  Horizontal Sampling Factor: " << (int) pCompo->horizontalSamplingFactor() << std::endl;
                    std::cerr << "  Vertical Sampling Factor: " << (int) pCompo->verticalSamplingFactor() << std::endl;
                    std::cerr << "  Quantization Table Destination Selector: " << (int) pCompo->quantizationTableDestinationSelector << std::endl;
#endif
                    compoMap[pCompo->componentsIdentifier] = i;
                    compos[i].setQuantizationSelector(pCompo->quantizationTableDestinationSelector);
                    compos[i].setHFactor(pCompo->horizontalSamplingFactor());
                    compos[i].setVFactor(pCompo->verticalSamplingFactor());

                    maxHFactor = std::max<int>(maxHFactor, pCompo->horizontalSamplingFactor());
                    maxVFactor = std::max<int>(maxVFactor, pCompo->verticalSamplingFactor());
                }
                mcuCountX = (mcuCountX + maxHFactor - 1) / maxHFactor;
                mcuCountY = (mcuCountY + maxVFactor - 1) / maxVFactor;
                mcuCount = mcuCountX * mcuCountY;
#ifdef JPEG_DUMP
                std::cerr << "MCU Count: " << (int) mcuCount << std::endl;
#endif
                for (int i = 0; i < pSOF->numberOfComponents; i++) {
                    int xi = ceil(width * frameCompoTable[i].horizontalSamplingFactor() / maxHFactor);
                    int yi = ceil(height * frameCompoTable[i].verticalSamplingFactor() / maxVFactor);
                    compos[i].resize(xi, yi, isProgress);
                }

                pData += 2 + parseBigendian16(pSOF->length);
                break;
            }
            case 0xDBFF: {
#ifdef JPEG_DUMP
                std::cerr << "Define Quantization Table" << std::endl;
                std::cerr << "-------------------------" << std::endl;
#endif
                uint8_t* temp = pData + 4;
                uint8_t* end = pData + parseBigendian16(pHeader->length) + 2;

                while (temp < end) {
                    JPEGQuantizationTableSpecification* p = reinterpret_cast<JPEGQuantizationTableSpecification*>(temp);
#ifdef JPEG_DUMP
                    std::cerr << "  Quantization Table:" << std::endl;
                    std::cerr << "  Quantization Table Element Precision: " << (int) p->quantizationTableElementPrecision() << std::endl;
                    std::cerr << "  Quantization Table Destination Identifier: " << (int) p->quantizationTableDestinationIdentifier() << std::endl;
#endif
                    int id = p->quantizationTableDestinationIdentifier();
                    int precision = p->quantizationTableElementPrecision();
                    temp += sizeof(JPEGQuantizationTableSpecification);
                    for (int i = 0; i < 64; i++) {
                        int qid = zz2quad[i];
                        uint16_t val = 0;
                        if (precision == 0) {
                            val = temp[i];
                        } else {
                            val = (temp[i << 1] << 8) | (temp[i << 1 | 1]);
                        }
                        quantizations[id][qid >> 3][qid & 0x07] = val;
                    }
                    temp += 64 * (precision + 1);
#ifdef JPEG_DUMP
                    std::cerr << "  Table:" << std::endl;
                    for (int i = 0; i < 8; i++) {
                        std::cerr << "    ";
                        for (int j = 0; j < 8; j++) {
                            std::cerr << quantizations[id][i][j] << " \n"[j == 7];
                        }
                    }
#endif
                }

                pData += 2 + parseBigendian16(pHeader->length);
                break;
            }
            case 0xC4FF: {
#ifdef JPEG_DUMP
                std::cerr << "Define Huffman Table" << std::endl;
                std::cerr << "--------------------" << std::endl;
#endif
                uint8_t* temp = pData + 4;
                uint8_t* end = pData + parseBigendian16(pHeader->length) + 2;

                while (temp < end) {
                    JPEGHuffmanTableSpecification* p = reinterpret_cast<JPEGHuffmanTableSpecification*>(temp);
#ifdef JPEG_DUMP
                    std::cerr << "  Huffman Table:" << std::endl;
                    std::cerr << "  Huffman Table Class: " << (int) p->tableClass() << std::endl;
                    std::cerr << "  Huffman Table Destination Identifier: " << (int) p->huffmanTableDestinationIdentifier() << std::endl;
                    std::cerr << "  Number of Huffman Codes of Length k: ";
                    for (int i = 0; i < 16; i++) std::cerr << (int) p->numberOfHuffmanCodesOfLengthK[i] << " ";
                    std::cerr << std::endl;
#endif
                    int tclass = p->tableClass();
                    int id = p->huffmanTableDestinationIdentifier();
                    temp += sizeof(JPEGHuffmanTableSpecification);
                    size_t parsedLength = huffmans[tclass << 2 | id].buildHuffmanTableFromDHT(p->numberOfHuffmanCodesOfLengthK, temp);
                    temp += parsedLength;
                }

                pData += 2 + parseBigendian16(pHeader->length);
                break;
            }
            case 0xDDFF: {
#ifdef JPEG_DUMP
                std::cerr << "Define Restart Interval" << std::endl;
                std::cerr << "-----------------------" << std::endl;
#endif
                JPEGDRI* pDRI = reinterpret_cast<JPEGDRI*>(pData);
                restartInterval = parseBigendian16(pDRI->restartInterval);
#ifdef JPEG_DUMP
                std::cerr << "Restart Interval: " << restartInterval << std::endl;
#endif
                pData += 2 + parseBigendian16(pDRI->length);
                break;
            }
            case 0xDAFF: {
#ifdef JPEG_DUMP
                std::cerr << "Start of Scan" << std::endl;
                std::cerr << "-------------" << std::endl;
#endif
                JPEGSOS* pSOS = reinterpret_cast<JPEGSOS*>(pData);
#ifdef JPEG_DUMP
                std::cerr << "Number of Image Components in Scan: " << (int) pSOS->numberOfComponentsInScan << std::endl;
#endif
                scanCompoTable.clear();
                JPEGScanSpecification* pScanSpec = reinterpret_cast<JPEGScanSpecification*>(pData + sizeof(JPEGSOS));
                for (int i = 0; i < pSOS->numberOfComponentsInScan; i++, pScanSpec++) {
                    scanCompoTable.push_back(*pScanSpec);
#ifdef JPEG_DUMP
                    std::cerr << "  Scan Component Selector: " << (int) pScanSpec->scanComponentsSelector << std::endl;
                    std::cerr << "  DC Entropy Coding Table Destination Selector: " << (int) pScanSpec->DCEntropyCodingTableDestinationSelector() << std::endl;
                    std::cerr << "  AC Entropy Coding Table Destination Selector: " << (int) pScanSpec->ACEntropyCodingTableDestinationSelector() << std::endl;
#endif
                    int cid = compoMap[pScanSpec->scanComponentsSelector];
                    compos[cid].setDChuffmanSelector(pScanSpec->DCEntropyCodingTableDestinationSelector());
                    compos[cid].setAChuffmanSelector(pScanSpec->ACEntropyCodingTableDestinationSelector());
                }
                this->pSOSEnd = reinterpret_cast<JPEGSOSEnd*>(pScanSpec);
#ifdef JPEG_DUMP
                std::cerr << "Start of Spectral of Predictor Selection: " << (int) pSOSEnd->startOfSpectralOfPredictorSelection << std::endl;
                std::cerr << "End of Spectral of Predictor Selection: " << (int) pSOSEnd->endOfSpectralOfPredictorSelection << std::endl;
                std::cerr << "Successive Approximation Bit Position High: " << (int) pSOSEnd->successiveApproximationBitPositionHigh() << std::endl;
                std::cerr << "Successive Approximation Bit Position Low: " << (int) pSOSEnd->successiveApproximationBitPositionLow() << std::endl;
#endif
                uint8_t* pECS = pData + 2 + parseBigendian16(pSOS->length);
                size_t parsedLength = parseECS(pECS, pEnd);

                pData += 2 + parseBigendian16(pSOS->length) + parsedLength;
                break;
            }
            case 0xD0FF:
            case 0xD1FF:
            case 0xD2FF:
            case 0xD3FF:
            case 0xD4FF:
            case 0xD5FF:
            case 0xD6FF:
            case 0xD7FF: {
                uint8_t* pECS = pData + 2;
                mcuIndexTemp = mcuIndex;
                if (isProgress) {
                    for (int i = 0; i < componentsCount; i++) {
                        compos[i].mcuIndexTemp = compos[i].mcuIndex;
                    }
                }
                size_t parsedLength = parseECS(pECS, pEnd);

                pData += 2 + parsedLength;
                break;
            }
            case 0xE0FF: {
#ifdef JPEG_DUMP
                std::cerr << "Application 0" << std::endl;
#endif
                pData += 2 + parseBigendian16(pHeader->length);
                break;
            }
            default: {
#ifdef JPEG_DUMP
                std::cerr << "Unknown Block" << std::endl;
                std::cerr << "-------------" << std::endl;
                std::cerr << "Marker: 0x" << std::hex << pHeader->marker << std::dec << std::endl;
#endif
                if ((pHeader->marker & 0xFF) != 0xFF) {
                    std::cerr << "Not a JPEG file" << std::endl;
                    delete[] filedata;
                    return Image(0, 0);
                }
                pData += 2 + parseBigendian16(pHeader->length);
                break;
            }
        }
    }

endOfParse:
    reconstruct(img);
    delete[] filedata;

    return img;
}

void JPEGUtility::write(const Image& img, const std::string& filename) {
    std::ofstream fp(filename, std::ios::binary);

    // SOI
    fp.put(0xFF); fp.put(0xD8);

    { // APP0
        JPEGMarkerHeader header;
        header.marker = 0xE0FF;
        header.length = getBigendian16(16);
        fp.write(reinterpret_cast<char*>(&header), sizeof(JPEGMarkerHeader)); // header
        fp.write("JFIF", 4); fp.put(0x00); // identifier
        fp.put(0x01); fp.put(0x02); // version
        fp.put(0x00); //density unit
        fp.write(getBigendian16PChar(1000), 2); // x density
        fp.write(getBigendian16PChar(1000), 2); // y density
        fp.put(0x00); // x thumbnail
        fp.put(0x00); // y thumbnail
    }

    { // SOF
        JPEGSOF sof;
        sof.marker = 0xC0FF;
        sof.length = getBigendian16(8 + 3 * 3);
        sof.samplePrecision = 8;
        sof.numberOfLines = getBigendian16(img.getHeight());
        sof.numberOfSamplesPerLine = getBigendian16(img.getWidth());
        sof.numberOfComponents = 3;
        fp.write(reinterpret_cast<char*>(&sof), sizeof(JPEGSOF));
        for (int i = 0; i < 3; i++) {
            JPEGFrameSpecification spec;
            spec.componentsIdentifier = i + 1;
            spec.samplingFactor = (1 << 4) | 1;
            spec.quantizationTableDestinationSelector = i ? 1 : 0;
            fp.write(reinterpret_cast<char*>(&spec), sizeof(JPEGFrameSpecification));
        }
    }

    { // DQT
        JPEGMarkerHeader header;
        header.marker = 0xDBFF;
        header.length = getBigendian16(2 + 2 * 65);
        fp.write(reinterpret_cast<char*>(&header), sizeof(JPEGMarkerHeader));

        // QTable for Y
        //    1 1 1 1 1 1 1 1
        //    1 1 1 1 1 1 1 1
        //    1 1 1 1 1 1 1 2
        //    1 1 1 1 1 1 2 2
        //    1 1 1 1 1 2 2 3
        //    1 1 1 1 2 2 3 3
        //    1 1 1 2 2 3 3 3
        //    1 1 2 2 3 3 3 3
        // 43-1, 11-2, 10-3
        fp.put(0x00);
        quantizations[0].create(8, 8, 0);
        for (int i = 0; i < 43; i++) {
            fp.put(0x01);
            int id = zz2quad[i];
            quantizations[0][id >> 3][id & 0x07] = 1;
        }
        for (int i = 0; i < 11; i++) {
            fp.put(0x02);
            int id = zz2quad[i + 43];
            quantizations[0][id >> 3][id & 0x07] = 2;
        }
        for (int i = 0; i < 10; i++) {
            fp.put(0x03);
            int id = zz2quad[i + 43 + 11];
            quantizations[0][id >> 3][id & 0x07] = 3;
        }

        // QTable for Cb and Cr
        //    1 1 1 1 2 3 3 3
        //    1 1 1 2 3 3 3 3
        //    1 1 1 3 3 3 3 3
        //    1 2 3 3 3 3 3 3
        //    2 3 3 3 3 3 3 3
        //    3 3 3 3 3 3 3 3
        //    3 3 3 3 3 3 3 3
        //    3 3 3 3 3 3 3 3
        // 10-1, 5-2, 49-3
        fp.put(0x01);
        quantizations[1].create(8, 8, 0);
        for (int i = 0; i < 10; i++) {
            fp.put(0x01);
            int id = zz2quad[i];
            quantizations[1][id >> 3][id & 0x07] = 1;
        }
        for (int i = 0; i < 5; i++) {
            fp.put(0x02);
            int id = zz2quad[i + 10];
            quantizations[1][id >> 3][id & 0x07] = 2;
        }
        for (int i = 0; i < 49; i++) {
            fp.put(0x03);
            int id = zz2quad[i + 10 + 5];
            quantizations[1][id >> 3][id & 0x07] = 3;
        }
    }

    { // DHT
        mcuCountX = (img.getWidth() + 7) >> 3;
        mcuCountY = (img.getHeight() + 7) >> 3;
        for (int i = 0; i < 3; i++) {
            compos[i].resize(mcuCountX * 8, mcuCountY * 8, true);
        }
        compos[0].setAChuffmanSelector(0);
        compos[0].setDChuffmanSelector(0);
        compos[1].setAChuffmanSelector(1);
        compos[1].setDChuffmanSelector(1);
        compos[2].setAChuffmanSelector(1);
        compos[2].setDChuffmanSelector(1);
        for (int i = 0; i < mcuCountY; i++) {
            for (int j = 0; j < mcuCountX; j++) {
                for (int y = 0; y < 8; y++) {
                    for (int x = 0; x < 8; x++) {
                        uint8_t R = img.getR(i * 8 + y, j * 8 + x);
                        uint8_t G = img.getG(i * 8 + y, j * 8 + x);
                        uint8_t B = img.getB(i * 8 + y, j * 8 + x);
                        uint8_t Y, Cb, Cr;
                        RGB2YCbCr(R, G, B, Y, Cb, Cr);
                        compos[0].setBlockAt(j, i, y, x, Y);
                        compos[1].setBlockAt(j, i, y, x, Cb);
                        compos[2].setBlockAt(j, i, y, x, Cr);
                    }
                }
            }
        }
        for (int c = 0; c < 3; c++) {
            compos[c].DCTAndQuant(quantizations[c ? 1 : 0]);
        }
        int pred[3] = {0, 0, 0};
        for (int i = 0; i < mcuCountY; i++) {
            for (int j = 0; j < mcuCountX; j++) {
                for (int c = 0; c < 3; c++) {
                    int DC = compos[c].getBlockAt(j, i, 0, 0);
                    compos[c].addBlockAt(j, i, 0, 0, -pred[c]);
                    pred[c] = DC;
                }
            }
        }

        JPEGMarkerHeader header;
        header.marker = 0xC4FF;
        int length = 2 + 17 * 4;
        uint8_t num[8][16];
        std::vector<uint8_t> data[8];
        // Y-DC HTable (0, 0)
        huffmans[0].buildSimpleDCHuffmanTable(num[0], data[0], codeMap[0]);
        for (int i = 0; i < 16; i++) length += num[0][i];
        // Cb/Cr-DC HTable (0, 1)
        huffmans[1].buildSimpleDCHuffmanTable(num[1], data[1], codeMap[1]);
        for (int i = 0; i < 16; i++) length += num[1][i];
        // Y-AC HTable (1, 0)
        huffmans[4].buildSimpleACHuffmanTable(num[4], data[4], codeMap[4]);
        for (int i = 0; i < 16; i++) length += num[4][i];
        // Cb/Cr-AC HTable (1, 1)
        huffmans[5].buildSimpleACHuffmanTable(num[5], data[5], codeMap[5]);
        for (int i = 0; i < 16; i++) length += num[5][i];

        header.length = getBigendian16(length);
        fp.write(reinterpret_cast<char*>(&header), sizeof(JPEGMarkerHeader));
        fp.put(0x00);
        fp.write(reinterpret_cast<char*>(num[0]), 16);
        fp.write(reinterpret_cast<char*>(data[0].data()), data[0].size());
        fp.put(0x01);
        fp.write(reinterpret_cast<char*>(num[1]), 16);
        fp.write(reinterpret_cast<char*>(data[1].data()), data[1].size());
        fp.put(0x10);
        fp.write(reinterpret_cast<char*>(num[4]), 16);
        fp.write(reinterpret_cast<char*>(data[4].data()), data[4].size());
        fp.put(0x11);
        fp.write(reinterpret_cast<char*>(num[5]), 16);
        fp.write(reinterpret_cast<char*>(data[5].data()), data[5].size());
    }

    { // SOS
        JPEGSOS sos;
        sos.marker = 0xDAFF;
        sos.length = getBigendian16(6 + 2 * 3);
        sos.numberOfComponentsInScan = 3;
        fp.write(reinterpret_cast<char*>(&sos), sizeof(JPEGSOS));
        for (int i = 0; i < 3; i++) {
            JPEGScanSpecification spec;
            spec.scanComponentsSelector = i + 1;
            spec.entropyCodingTableDestinationSelector = i ? 0x11 : 0x00;
            fp.write(reinterpret_cast<char*>(&spec), sizeof(JPEGScanSpecification));
        }
        JPEGSOSEnd sosEnd;
        sosEnd.startOfSpectralOfPredictorSelection = 0;
        sosEnd.endOfSpectralOfPredictorSelection = 63;
        sosEnd.successiveApproximationBitPosition = 0;
        fp.write(reinterpret_cast<char*>(&sosEnd), sizeof(JPEGSOSEnd));
    }

    { // ECS
        scanData.clear();
        size_t byteOffset = 0, bitOffset = 0;
        for (int i = 0; i < mcuCountY; i++) {
            for (int j = 0; j < mcuCountX; j++) {
                for (int c = 0; c < 3; c++) {
                    { // DC
                        int rs = getCodeLength(compos[c].getBlockAt(j, i, 0, 0));
                        auto code = codeMap[compos[c].DCHuffman()][rs];
                        pushValue(byteOffset, bitOffset, code);
                        int ext = getExtendCode(compos[c].getBlockAt(j, i, 0, 0));
                        pushValue(byteOffset, bitOffset, std::make_pair(rs, ext));
                    }
                    // AC
                    int list[64];
                    for (int y = 0; y < 8; y++) {
                        for (int x = 0; x < 8; x++) {
                            list[quad2zz[y * 8 + x]] = compos[c].getBlockAt(j, i, y, x);
                        }
                    }
                    int last = 63, zeroCnt = 0;
                    for (; last >= 1; last--) if (list[last]) break;
                    for (int aid = 1; aid <= last; aid++) {
                        if (list[aid]) {
                            int s = getCodeLength(list[aid]);
                            int rs = (zeroCnt << 4) | s;
                            auto code = codeMap[compos[c].ACHuffman()][rs];
                            pushValue(byteOffset, bitOffset, code);
                            int ext = getExtendCode(list[aid]);
                            pushValue(byteOffset, bitOffset, std::make_pair(s, ext));
                            zeroCnt = 0;
                        } else {
                            ++zeroCnt;
                            if (zeroCnt >= 16) { // ZRL
                                auto code = codeMap[compos[c].ACHuffman()][0xF0];
                                pushValue(byteOffset, bitOffset, code);
                                zeroCnt -= 16;
                            }
                        }
                    }
                    if (last < 63) { // EOB
                        auto code = codeMap[compos[c].ACHuffman()][0x00];
                        pushValue(byteOffset, bitOffset, code);
                    }
                }
            }
        }
        if (bitOffset) {
            scanData[byteOffset] |= ((1 << (8 - bitOffset)) - 1);
            ++byteOffset;
            bitOffset = 0;
        } else {
            scanData.pop_back();
        }
        assert(scanData.size() == byteOffset);
        std::vector<uint8_t> dataWithBitstuff;
        for (int i = 0, j = 0; i < byteOffset; i++) {
            dataWithBitstuff.push_back(scanData[i]);
            if (scanData[i] == 0xFF) dataWithBitstuff.push_back(0x00);
        }
        fp.write(reinterpret_cast<char*>(dataWithBitstuff.data()), dataWithBitstuff.size());
    }

    // EOI
    fp.put(0xFF); fp.put(0xD9);

    fp.close();
}
