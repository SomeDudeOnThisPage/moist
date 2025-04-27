#include "tiff_data.hpp"

#include <format>
#include <iostream>
#include <string>

// TODO: some variant crap to enable 8,32 bit etc... test-data for now uses 16 bit samples.
incremental_meshing::TiffData::TiffData(const std::string& pattern, const uint32_t first_file, const uint32_t num_files) : data(std::vector<uint16_t>()), _depth(num_files)
{
    // unfortunately, often tiff data isn't set up as a directory, so we need to iterate and open each file manually to be sure...
    data.resize(_depth);
    for (auto n = 0; n < num_files; n++)
    {
        try
        {
            const std::string filename = std::vformat(pattern, std::make_format_args(n + first_file));
            auto tiff = std::unique_ptr<TIFF, decltype(&TIFFClose)>(
                TIFFOpen(filename.c_str(), "r"),
                &TIFFClose
            );

            uint32_t sample_bits = 16;
            TIFFGetField(tiff.get(), TIFFTAG_IMAGEWIDTH, &_width);
            TIFFGetField(tiff.get(), TIFFTAG_IMAGELENGTH, &_height);
            TIFFGetField(tiff.get(), TIFFTAG_BITSPERSAMPLE, &sample_bits);

            if (n == 0)
            {
                data.resize(_width * _height * _depth);
            }

            OOC_DEBUG("loading .tiff: w = " << _width << ", h = " << _height << ", sample_bits = " << sample_bits);

            //_data[n].resize(_height);
            for (auto row = 0; row < _height; row++)
            {
                //_data[n][row].resize(_width);
                if (TIFFReadScanline(tiff.get(), &data[n * _width *_height + row * _width] /*_data[n][row].data()*/, row, 0) < 0) {
                    OOC_ERROR("Failed to read scanline " << row);
                }
            }
        }
        catch (const std::format_error& e)
        {
            OOC_ERROR("invalid file pattern: " << e.what());
        }
    }
}
