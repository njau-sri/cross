#include <iomanip>
#include <iostream>
#include "cmdline.h"

CmdLine::CmdLine(const std::string &info)
    : info_(info)
{
}

void CmdLine::add(const std::string &arg, const std::string &info)
{
    flags_[arg] = false;
    flags_info_[arg] = info;
}

void CmdLine::add(const std::string &arg, const std::string &info, const std::string &val)
{
    args_[arg] = val;
    args_info_[arg] = info;
}

void CmdLine::help() const
{
    std::cerr << std::left;
    std::cerr << "usage: " << info_ << "\n";

    for (auto &e : args_)
        std::cerr << "  " << std::setw(25) << e.first + "  <" + e.second + ">" << " " << args_info_.at(e.first) << "\n";

    for (auto &e : flags_)
        std::cerr << "  " << std::setw(25) << e.first << " " << flags_info_.at(e.first) << "\n";
}

void CmdLine::parse(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i) {
        if (args_.count(argv[i]) != 0) {
            if (++i != argc && args_.count(argv[i]) == 0 && flags_.count(argv[i]) == 0)
                args_[argv[i - 1]] = argv[i];
            else
                std::cerr << "ERROR: missing an argument for: " << argv[--i] << "\n";
        }
        else if (flags_.count(argv[i]) != 0) {
            flags_[argv[i]] = true;
        }
        else {
            std::cerr << "ERROR: unrecognized command line argument: " << argv[i] << "\n";
        }
    }
}

bool CmdLine::has(const std::string &arg) const
{
    return flags_.at(arg);
}

std::string CmdLine::get(const std::string &arg) const
{
    return args_.at(arg);
}
