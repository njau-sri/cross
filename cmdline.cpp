#include <iomanip>
#include <iostream>
#include "cmdline.h"

CmdLine::CmdLine(const std::string &usage)
    : m_usage(usage)
{
}

void CmdLine::add(const std::string &arg, const std::string &info)
{
    m_flagval[arg] = false;
    m_flaginfo[arg] = info;
}

void CmdLine::add(const std::string &arg, const std::string &info, const std::string &val)
{
    m_argval[arg] = val;
    m_arginfo[arg] = info;
}

void CmdLine::help() const
{
    std::cerr << std::left;
    std::cerr << "usage: " << m_usage << "\n";

    for (auto &e : m_argval)
        std::cerr << "  " << std::setw(25) << e.first + "  <" + e.second + "> " << m_arginfo.at(e.first) << "\n";

    for (auto &e : m_flagval)
        std::cerr << "  " << std::setw(25) << e.first << " " << m_flaginfo.at(e.first) << "\n";
}

void CmdLine::parse(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i) {
        if (m_argval.count(argv[i]) != 0) {
            if (++i != argc && m_argval.count(argv[i]) == 0 && m_flagval.count(argv[i]) == 0)
                m_argval[argv[i - 1]] = argv[i];
            else
                std::cerr << "ERROR: missing an argument for: " << argv[--i] << "\n";
        }
        else if (m_flagval.count(argv[i]) != 0) {
            m_flagval[argv[i]] = true;
        }
        else {
            std::cerr << "ERROR: unrecognized command line argument: " << argv[i] << "\n";
        }
    }
}

bool CmdLine::has(const std::string &arg) const
{
    return m_flagval.at(arg);
}

std::string CmdLine::get(const std::string &arg) const
{
    return m_argval.at(arg);
}
