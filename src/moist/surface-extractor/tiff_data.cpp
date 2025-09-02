#include "tiff_data.hpp"

#include <format>
#include <iostream>
#include <string>

#include "moist/core/timer.hpp"

moist::Tiff::Tiff(const std::string& pattern, const uint32_t first_file, const uint32_t num_files) : data(std::vector<float>()), _depth(num_files)
{
    moist::ScopeTimer timer("Tiff::Tiff");
    for (auto n = 0; n < num_files; n++)
    {
        try
        {
            const std::string filename = std::vformat(pattern, std::make_format_args(n + first_file));
            auto tiff = std::unique_ptr<TIFF, decltype(&TIFFClose)>( // wrap in unique-pointer for automatic closing...
                TIFFOpen(filename.c_str(), "r"),
                &TIFFClose
            );

            TIFFGetField(tiff.get(), TIFFTAG_IMAGEWIDTH, &_width);
            TIFFGetField(tiff.get(), TIFFTAG_IMAGELENGTH, &_height);
            TIFFGetField(tiff.get(), TIFFTAG_BITSPERSAMPLE, &_sample_bits);

            if (n == 0)
            {
                this->data.reserve(this->_width * this->_height * this->_depth);
            }

            OOC_DEBUG("loading " << filename << ": w = " << _width << ", h = " << _height << ", sample_bits = " << _sample_bits);

            _vector line_data;
            float max;
            switch (_sample_bits)
            {
                case 8:
                    line_data = std::vector<uint8_t>(this->_width);
                    max = static_cast<float>(std::numeric_limits<uint8_t>::max());
                    break;
                case 16:
                    line_data = std::vector<uint16_t>(this->_width);
                    max = static_cast<float>(std::numeric_limits<uint16_t>::max());
                    break;
                case 32:
                    line_data = std::vector<uint32_t>(this->_width);
                    max = static_cast<float>(std::numeric_limits<uint32_t>::max());
                    break;
                case 64:
                    line_data = std::vector<uint64_t>(this->_width);
                    max = static_cast<float>(std::numeric_limits<uint64_t>::max());
                    break;
                default:
                    OOC_ERROR("Invalid tiff sample bits: " << _sample_bits);
            }

            for (auto row = 0; row < _height; row++)
            {
                std::visit([this, &tiff, &n, &row, &max](auto& vector)
                {
                    // using T = typename std::decay_t<decltype(vector)>::value_type;

                    if (TIFFReadScanline(tiff.get(), vector.data(), row, 0) < 0)
                    {
                        OOC_ERROR("Failed to read scanline " << row);
                    }

                    // enable parallel insertion into the vector by not using push_back
                    for (uint32_t i = 0; i < vector.size(); i++)
                    {
                        this->data[n * _width *_height + row * _width + i] = static_cast<float>(vector[i]) > 0 ? 1.0f : 0.0f;
                    }
                }, line_data);
            }
        }
        catch (const std::format_error& e)
        {
            OOC_ERROR("invalid file pattern: " << e.what());
        }
    }
}
