#ifndef __TIFF_DATA_HPP
#define __TIFF_DATA_HPP

#include <variant>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <iostream>

#include <tiffio.h>

#include "../core.hpp"

namespace ooc
{
    class TiffData
    {
    public:
        TiffData(std::filesystem::path first_file, uint32_t n_files);
        ~TiffData() = default;

        uint32_t width() { return _width; };
        uint32_t height() { return _height; };
        uint32_t depth() { return _depth; };

        std::vector<std::vector<std::vector<uint16_t>>>& data() { return _data; };
        std::vector<std::vector<std::vector<uint16_t>>> _data;
    private:
        uint32_t _width;
        uint32_t _height;
        uint32_t _depth;
    };
}

#endif // __TIFF_DATA_HPP
