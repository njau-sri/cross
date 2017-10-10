#pragma once

#include <map>
#include <string>

class CmdLine
{
public:
    CmdLine(const std::string &info);

    void add(const std::string &arg, const std::string &info);

    void add(const std::string &arg, const std::string &info, const std::string &val);

    void help() const;

    void parse(int argc, char *argv[]);

    bool has(const std::string &arg) const;

    std::string get(const std::string &arg) const;

private:
    std::map<std::string, std::string> args_;
    std::map<std::string, std::string> args_info_;
    std::map<std::string, bool> flags_;
    std::map<std::string, std::string> flags_info_;
    std::string info_;
};
