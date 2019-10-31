#pragma once

#include "Image.hpp"
#include <cstdint>
#include <string>
#include "Mat.hpp"

struct HuffmanNode {
    std::shared_ptr<HuffmanNode> c[2];
    uint16_t value;
    bool isLeaf = false;

    HuffmanNode() = default;
    HuffmanNode(uint16_t value) : c(), value(value), isLeaf(true) {}
};

class HuffmanTable {
public:
    HuffmanTable() { root = std::make_shared<HuffmanNode>(); }
    ~HuffmanTable() = default;

    size_t buildHuffmanTableFromDHT(const uint8_t* num, const uint8_t* data);
    void buildSimpleDCHuffmanTable(uint8_t* num, std::vector<uint8_t>& data, std::pair<int, int>* codeMap);
    void buildSimpleACHuffmanTable(uint8_t* num, std::vector<uint8_t>& data, std::pair<int, int>* codeMap);
    uint16_t parse(const uint8_t* data, size_t size, size_t& byteOffset, size_t& bitOffset) const;
    void getCodeMap(std::pair<int, int>* codeMap);

private:
    std::shared_ptr<HuffmanNode> root;
};

#pragma pack(push, 1)

struct JPEGSOI {
    uint16_t SOI;
};

struct JPEGMarkerHeader {
    uint16_t marker;
    uint16_t length;
};

struct JPEGSOF : JPEGMarkerHeader {
    uint8_t samplePrecision;
    uint16_t numberOfLines;
    uint16_t numberOfSamplesPerLine;
    uint8_t numberOfComponents;
};

struct JPEGFrameSpecification {
    uint8_t componentsIdentifier;
    uint8_t samplingFactor;
    uint8_t quantizationTableDestinationSelector;
    uint8_t horizontalSamplingFactor() const { return samplingFactor >> 4u; }
    uint8_t verticalSamplingFactor() const { return samplingFactor & 0x0Fu; }
};

struct JPEGSOS : JPEGMarkerHeader {
    uint8_t numberOfComponentsInScan;
};

struct JPEGScanSpecification {
    uint8_t scanComponentsSelector;
    uint8_t entropyCodingTableDestinationSelector;
    uint8_t DCEntropyCodingTableDestinationSelector() const { return entropyCodingTableDestinationSelector >> 4u; }
    uint8_t ACEntropyCodingTableDestinationSelector() const { return entropyCodingTableDestinationSelector & 0x0Fu; }
};

struct JPEGSOSEnd {
    uint8_t startOfSpectralOfPredictorSelection;
    uint8_t endOfSpectralOfPredictorSelection;
    uint8_t successiveApproximationBitPosition;
    uint8_t successiveApproximationBitPositionHigh() const { return successiveApproximationBitPosition >> 4u; }
    uint8_t successiveApproximationBitPositionLow() const { return successiveApproximationBitPosition & 0x0Fu; }
};

struct JPEGQuantizationTableSpecification {
    uint8_t data;
    uint8_t quantizationTableElementPrecision() const { return data >> 4u; }
    uint8_t quantizationTableDestinationIdentifier() const { return data & 0x0Fu; }
};

struct JPEGHuffmanTableSpecification {
    uint8_t data;
    uint8_t numberOfHuffmanCodesOfLengthK[16];
    uint8_t tableClass() const { return data >> 4u; }
    uint8_t huffmanTableDestinationIdentifier() const { return data & 0x0Fu; }
};

struct JPEGDRI : JPEGMarkerHeader {
    uint16_t restartInterval;
};

#pragma pack(pop)

class JPEGUtility {
public:
    static JPEGUtility* JPEG() {
        if (!singleton) singleton = new JPEGUtility();
        return singleton;
    }

    Image read(const std::string& filename);
    void write(const Image& img, const std::string& filename);

private:
    JPEGUtility() {
        for (int i = 0; i < 4; i++) quantizations[i].create(8, 8, 0);
        std::fill(compoMap, compoMap + 256, 0);
    }
    ~JPEGUtility() = default;

    static JPEGUtility* singleton;

    Mat<uint16_t> quantizations[4];
    HuffmanTable huffmans[8];
    std::pair<int, int> codeMap[8][256];

    size_t width;
    size_t height;
    size_t componentsCount;
    size_t mcuIndex;
    size_t mcuIndexTemp;
    size_t mcuCount;
    size_t mcuCountX;
    size_t mcuCountY;
    size_t restartInterval = 0;
    uint8_t samplePrecision;
    bool isProgress;

    std::vector<JPEGFrameSpecification> frameCompoTable;
    std::vector<JPEGScanSpecification> scanCompoTable;
    int compoMap[256];
    JPEGSOSEnd* pSOSEnd;

    std::vector<uint8_t> scanData;
    int predDC[4];

    constexpr static int quad2zz[] = {
            0,  1,  5,  6, 14, 15, 27, 28,
            2,  4,  7, 13, 16, 26, 29, 42,
            3,  8, 12, 17, 25, 30, 41, 43,
            9, 11, 18, 24, 31, 40, 44, 53,
            10, 19, 23, 32, 39, 45, 52, 54,
            20, 22, 33, 38, 46, 51, 55, 60,
            21, 34, 37, 47, 50, 56, 59, 61,
            35, 36, 48, 49, 57, 58, 62, 63
    };
    constexpr static int zz2quad[] = {
            0,  1,  8, 16,  9,  2,  3, 10,
            17, 24, 32, 25, 18, 11,  4,  5,
            12, 19, 26, 33, 40, 48, 41, 34,
            27, 20, 13,  6,  7, 14, 21, 28,
            35, 42, 49, 56, 57, 50, 43, 36,
            29, 22, 15, 23, 30, 37, 44, 51,
            58, 59, 52, 45, 38, 31, 39, 46,
            53, 60, 61, 54, 47, 55, 62, 63
    };

    int getExtendValue(size_t& byteOffset, size_t& bitOffset, size_t length = 1);
    void pushValue(size_t& byteOffset, size_t& bitOffset, const std::pair<int, int>& code);
    int getCodeLength(int val);
    int getExtendCode(int val);
    size_t parseECS(uint8_t* data, uint8_t* end);
    void parseECS_base();
    Mat<int> parseECSBlock_base(int cid, size_t& byteOffset, size_t& bitOffset);
    void parseECS_prog_DC();
    void parseECS_prog_AC();
    void reconstruct(Image& img);
};
