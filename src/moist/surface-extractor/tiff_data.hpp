#ifndef MOIST_SURFACE_EXTRACTOR_TIFF_HPP_
#define MOIST_SURFACE_EXTRACTOR_TIFF_HPP_

#include <variant>
#include <vector>
#include <cstdint>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <iostream>

#include <tiffio.h>

#include "moist/core/defines.hpp"

namespace moist
{
    using _vector = std::variant<
        std::vector<uint8_t>,
        std::vector<uint16_t>,
        std::vector<uint32_t>,
        std::vector<uint64_t>
    >;

    class Tiff
    {
    public:
        Tiff(const std::string& pattern, uint32_t first_file, uint32_t num_files);
        ~Tiff() = default;

        uint32_t width() const { return _width; };
        uint32_t height() const { return _height; };
        uint32_t depth() const { return _depth; };

        uint8_t bits() const { return _sample_bits; }

        std::vector<float> data;
    private:
        uint32_t _width;
        uint32_t _height;
        uint32_t _depth;
        uint8_t _sample_bits;
    };
}

#endif // MOIST_SURFACE_EXTRACTOR_TIFF_HPP_
