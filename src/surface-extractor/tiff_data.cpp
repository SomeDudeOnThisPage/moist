#include "tiff_data.hpp"

incremental_meshing::TiffData::TiffData(std::filesystem::path first_file, uint32_t n_files) : _depth(n_files)
{
    // TODO: some variant crap to enable 8,32 bit etc... test-data for now uses 16 bit samples.
    _data = std::vector<std::vector<std::vector<uint16_t>>>();

    auto tiff = std::unique_ptr<TIFF, decltype(&TIFFClose)>(
        TIFFOpen(first_file.c_str(), "r"),
        &TIFFClose
    );

    uint32_t sample_bits = 16;
    TIFFGetField(tiff.get(), TIFFTAG_IMAGEWIDTH, &_width);
    TIFFGetField(tiff.get(), TIFFTAG_IMAGELENGTH, &_height);
    TIFFGetField(tiff.get(), TIFFTAG_BITSPERSAMPLE, &sample_bits);

    OOC_DEBUG("loaded .tiff: w = " << _width << ", h = " << _height << ", sample_bits = " << sample_bits);

    auto directory = std::vector<std::vector<uint16_t>>();
    directory.resize(_width);
    std::cout << directory.size() << std::endl;
    for (auto row = 0; row < _height; row++)
    {
        directory[row].resize(_width);
        if (TIFFReadScanline(tiff.get(), directory[row].data(), row, 0) < 0) {
            OOC_ERROR("Failed to read scanline " << row);
        }
    }

    this->_data.resize(_depth);
    this->_data[0] = directory;
}
