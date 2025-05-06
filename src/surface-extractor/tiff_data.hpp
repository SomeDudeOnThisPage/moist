#ifndef __TIFF_DATA_HPP
#define __TIFF_DATA_HPP

#include <variant>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <iostream>

#include <tiffio.h>

#include "core.hpp"

namespace incremental_meshing
{
    class TiffData
    {
    public:
        TiffData(const std::string& pattern, uint32_t first_file, uint32_t num_files);
        ~TiffData() = default;

        uint32_t width() const { return _width; };
        uint32_t height() const { return _height; };
        uint32_t depth() const { return _depth; };

        uint8_t bits() const { return _sample_bits; }

        std::vector<uint16_t> data;
    private:
        uint32_t _width;
        uint32_t _height;
        uint32_t _depth;
        uint8_t _sample_bits;
    };
}

#endif // __TIFF_DATA_HPP
