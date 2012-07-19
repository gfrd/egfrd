#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <string>
#include <exception>
#include <stdexcept>

class illegal_state: public std::exception
{
public:
    illegal_state(std::string const& str): str_(str) {}

    virtual ~illegal_state() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class illegal_argument: public std::exception
{
public:
    illegal_argument(std::string const& str): str_(str) {}

    virtual ~illegal_argument() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class not_found: public std::exception
{
public:
    not_found(std::string const& str): str_(str) {}

    virtual ~not_found() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class already_exists: public std::exception
{
public:
    already_exists(std::string const& str): str_(str) {}

    virtual ~already_exists() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class unsupported: public std::exception
{
public:
    unsupported(std::string const& str): str_(str) {}

    virtual ~unsupported() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class propagation_error: public std::runtime_error
{
public:
    propagation_error(std::string const& msg): std::runtime_error(msg) {}

    virtual ~propagation_error() throw() {}

private:
    std::string str_;
};

class illegal_propagation_attempt: public std::runtime_error
{
public:
    illegal_propagation_attempt(std::string const& msg): std::runtime_error(msg) {}
    // TODO: remove string output once everything works

    virtual ~illegal_propagation_attempt() throw() {}

private:
    std::string str_;
};

class not_implemented: public std::exception
{
public:
    not_implemented(std::string const& str): str_(str) {}

    virtual ~not_implemented() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class no_space: public std::exception
{
public:
    no_space(std::string const& str = ""): str_(str) {}

    virtual ~no_space() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

#endif /* EXCEPTIONS_HPP */
