#include "../include/vector_exception.h"

vector_exception::vector_exception(const std::string& msg) 
    : message(msg) {}

const char* vector_exception::what() const noexcept{
    return message.c_str();
}