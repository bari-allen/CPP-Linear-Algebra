#ifndef VECTOR_EXCEPTION_H
#define VECTOR_EXCEPTION_H

#include <exception>
#include <string>

class vector_exception : public std::exception {
    public:
        explicit vector_exception(const std::string& msg);

        const char* what() const noexcept override;

    private:
        std::string message;
};

#endif