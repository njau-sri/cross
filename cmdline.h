#pragma once

#include <map>
#include <string>

class CmdLine
{
public:
    CmdLine(const std::string &usage);

    void add(const std::string &arg, const std::string &info);

    void add(const std::string &arg, const std::string &info, const std::string &val);

    void help() const;

    void parse(int argc, char *argv[]);

    bool has(const std::string &arg) const;

    std::string get(const std::string &name) const;

private:
    std::map<std::string, std::string> m_argval;
    std::map<std::string, std::string> m_arginfo;
    std::map<std::string, bool> m_flagval;
    std::map<std::string, std::string> m_flaginfo;
    std::string m_usage;
};
