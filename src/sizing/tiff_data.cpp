#include "tiff_data.hpp"

ooc::TiffData::TiffData(std::filesystem::path first_file)
{
    this->_data = std::vector<std::vector<std::vector<uint16_t>>>();

    auto tiff = std::unique_ptr<TIFF, decltype(&TIFFClose)>(
        TIFFOpen(first_file.c_str(), "r"),
        &TIFFClose
    );

    uint32_t sample_bits;
    TIFFGetField(tiff.get(), TIFFTAG_IMAGEWIDTH, &this->_width);
    TIFFGetField(tiff.get(), TIFFTAG_IMAGELENGTH, &this->_height);
    TIFFGetField(tiff.get(), TIFFTAG_BITSPERSAMPLE, &sample_bits);

    OOC_DEBUG("tiff w = " << this->_width << ", h = " << this->_height << ", sample_bits = " << sample_bits);

    auto directory = std::vector<std::vector<uint16_t>>();
    directory.resize(this->_width);
    std::cout << directory.size() << std::endl;
    for (auto row = 0; row < this->_height; row++)
    {
        directory[row].resize(this->_width);
        if (TIFFReadScanline(tiff.get(), directory[row].data(), row, 0) < 0) {
            OOC_ERROR("Failed to read scanline " << row);
        }
    }

    //for (int i = 0; i < _height; i++)
    //{
    //    for (int j = 0; j < _width; j++)
    //    {
    //        std::cout << i << "," << j << " = " << directory[j][i] << std::endl;
    //    }
    //}
    this->_data.resize(1);
    this->_data[0] = directory;
}
