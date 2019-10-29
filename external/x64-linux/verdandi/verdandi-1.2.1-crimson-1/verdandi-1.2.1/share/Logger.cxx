// Copyright (C) 2008, INRIA
// Author(s): Marc Fragu
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_SHARE_LOGGER_CXX

#include "Logger.hxx"


namespace Verdandi
{


    ///////////////////
    // STATIC FIELDS //
    ///////////////////

    bool Logger::is_initialized_ = false;
    bool Logger::is_active_ = VERDANDI_LOG_IS_ACTIVE;
    string Logger::file_name_;
    int Logger::options_ = 0;
    int Logger::default_options_ = 0;
    int Logger::logging_level_;
    Logger::CommandMap Logger::command_;


    ////////////////////
    // STATIC METHODS //
    ////////////////////


    //! Initializes the Logger class.
    void Logger::Initialize()
    {
        InitializeDefaultOptions();
        InitializeOptions();
        InitializeFilename();
        InitializeLevel();
        InitializeCommand();
        EmptyFile();
        is_initialized_ = true;
    }



    //! Initializes the Logger class.
    /*!
      \param[in] configuration_file the configuration file.
      \param[in] section_name the section in \a configuration_file where the
      configuration is to be read.
    */
    void Logger::Initialize(string configuration_file, string section_name)
    {
        InitializeDefaultOptions(configuration_file, section_name);
        InitializeOptions();
        InitializeFilename(configuration_file, section_name);
        InitializeLevel(configuration_file, section_name);
        InitializeCommand();
        EmptyFile();
        is_initialized_ = true;
    }


    //! Clears the Logger class.
    void Logger::Finalize()
    {
        file_name_ = "";
        options_ = 0;
        default_options_ = 0;
        command_.clear();
    }


    //! Sets 'options_' to 'default_options_'.
    void Logger::InitializeOptions()
    {
        options_ = default_options_;
    }


    //! Sets the name of the log file.
    /*!
      \param[in] file_name path to the log file.
    */
    void Logger::SetFileName(string file_name)
    {
        file_name_ = file_name;
#if defined(VERDANDI_WITH_MPI)
        file_name_ = find_replace(file_name_, "%{rank}",
                                  to_str(MPI::COMM_WORLD.Get_rank()));
#else
        file_name_ = find_replace(file_name_, "%{rank}", "");
#endif
    }


    //! Activates or deactivates a specific option.
    /*!
      \param[in] option option to (des)activate.
      \param[in] value value.
    */
    void Logger::SetOption(int option, bool value)
    {
        if (!CheckStatus())
            return;
        if (value)
            options_ = options_ | option;
        else
            options_ = options_ & (~option);
    }


    //! Activates or deactivates "stdout" option.
    /*!
      \param[in] value value.
    */
    void Logger::SetStdout(bool value)
    {
        if (!CheckStatus())
            return;
        SetOption(stdout_, value);
    }


    //! Activates or deactivates "file" option.
    /*!
      \param[in] value value.
    */
    void Logger::SetFile(bool value)
    {
        if (!CheckStatus())
            return;
        SetOption(file_, value);
    }


    //! Activates or deactivates "uppercase" option.
    /*!
      \param[in] value value.
    */
    void Logger::SetUppercase(bool value)
    {
        if (!CheckStatus())
            return;
        SetOption(uppercase_, value);
    }


    //! Sets the level of verbosity.
    /*!
      \param[in] level level of verbosity.
    */
    void Logger::SetLoggingLevel(int level)
    {
        if (!CheckStatus())
            return;
        logging_level_ = level;
    }


    //! Activates the Logger.
    void Logger::Activate()
    {
        is_active_ = true;
    }


    //! Deactivates the Logger.
    void Logger::Deactivate()
    {
        is_active_ = false;
    }


    //! Writes a message in the log file.
    /*!
      \tparam S type of the message, which must be convertible to a string
      through 'ostringstream& operator << (ostringstream&, S& message)'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <int LEVEL, class T, class S>
    void Logger::Log(const T& object, const S& message, int options)
    {
        if (!CheckStatus(options))
            return;
        if (LEVEL >= logging_level_)
            Log(object, to_str(message), options);
    }


    //! Writes a message in the log file.
    /*!
      \tparam S type of the message, which must be convertible to a string
      through 'ostringstream& operator << (ostringstream&, S& message)'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <class T, class S>
    void Logger::Log(const T& object, const S& message, int options)
    {
        Log(object, to_str(message), options);
    }


    //! Writes a message in the log file.
    /*!
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <class T>
    void Logger::Log(const T& object, string message, int options)
    {
        if (!CheckStatus(options))
            return;
        WriteMessage(object, message, options);
    }


    //! Writes a message in the standard output and in the log file.
    /*! The message is always sent to the standard output, and it is possibly
      written in a log file if the logging level is lower than or equal to
      'VERDANDI_STDOUT_LOGGING_LEVEL'.
      \tparam S type of the message, which must be convertible to a string
      through 'ostringstream& operator << (ostringstream&, S& message)'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
    */
    template <class T, class S>
    void Logger::StdOut(const T& object, const S& message)
    {
        StdOut(object, to_str(message));
    }


    //! Writes a message in the standard output and in the log file.
    /*! The message is always sent to the standard output, and it is possibly
      written in a log file if the logging level is lower than or equal to
      'VERDANDI_STDOUT_LOGGING_LEVEL'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
    */
    template <class T>
    void Logger::StdOut(const T& object, string message)
    {
        WriteMessage(object, message, stdout_);
        if (VERDANDI_STDOUT_LOGGING_LEVEL >= logging_level_)
            WriteMessage(object, message, options_ & ~stdout_);
    }


    //! Reads a specific command and applies the corresponding treatment.
    /*!
      \param[in] command the specific command.
      \param[in] parameter the parameter of the command.
      \param[in] options options.
      \note Only 'hline' command is yet supported.
    */
    void Logger::Command(string command,
                         string parameter, int options)
    {
        if (!CheckStatus(options))
            return;
        CommandMap::const_iterator im;
        im = command_.find(command);
        LogCommand function_pointer;
        Logger logger;
        if (im != command_.end())
        {
            function_pointer = im->second;
            (logger.*function_pointer)(parameter, options);
        }
    }


    //! Reads a specific command and applies the corresponding treatment.
    /*! The message is always sent to the standard output, and it is possibly
      written in a log file if the logging level is lower than or equal to
      'VERDANDI_STDOUT_LOGGING_LEVEL'.
      \param[in] command the specific command.
      \param[in] parameter the parameter of the command.
      \note Only 'hline' command is yet supported.
    */
    void Logger::StdOutCommand(string command, string parameter)
    {
        Command(command, parameter, stdout_);
        if (VERDANDI_STDOUT_LOGGING_LEVEL >= logging_level_)
            Command(command, parameter, options_ & ~stdout_);
    }


    ////////////////////////////
    // PRIVATE STATIC METHODS //
    ////////////////////////////


    //! Initializes the default options.
    void Logger::InitializeDefaultOptions()
    {
        default_options_ =  VERDANDI_LOG_OPTIONS;
    }


    //! Initializes the default options.
    /*!
      \param[in] configuration_file path to the configuration file.
      \param[in] section_name the section in \a configuration_file where the
      configuration is to be read.
    */
    void Logger::InitializeDefaultOptions(string configuration_file,
                                          string section_name)
    {
        // To store the existing options and their corresponding "bit
        // position" in 'options_'.
        map<string, int> valid_options;
        valid_options["stdout"] = stdout_;
        valid_options["file"] = file_;
        valid_options["uppercase"] = uppercase_;

        VerdandiOps configuration(configuration_file);
        configuration.SetPrefix(section_name);
        string tmp;
        configuration.Set("default_options", "", "", tmp);

        vector<string> options_vector = split(tmp);
        if (options_vector.size() == 0)
        {
            InitializeDefaultOptions();
            return;
        }
        map<string, int>::const_iterator im;
        for (unsigned int i = 0; i < options_vector.size(); i++)
        {
            im = valid_options.find(options_vector[i]);
            if (im == valid_options.end())
            {
                InitializeDefaultOptions();
                return;
            }
            default_options_ = default_options_ | im->second;
        }
    }


    //! Initializes the name of the log file.
    void Logger::InitializeFilename()
    {
        // If the file name is not empty, it has already been set using
        // 'SetFileName' and it should not be overwritten.
        if (file_name_.empty())
        {
            file_name_ = find_replace(VERDANDI_LOG_FILENAME, "%{D}",
                                      GenerateDate());
#if defined(VERDANDI_WITH_MPI)
            file_name_ = find_replace(file_name_, "%{rank}",
                                      to_str(MPI::COMM_WORLD.Get_rank()));
#else
            file_name_ = find_replace(file_name_, "%{rank}", "");
#endif
        }
    }


    //!  Initializes the name of the the log file.
    /*!
      \param[in] configuration_file path to the configuration file.
      \param[in] section_name the section in \a configuration_file where the
      configuration is to be read.
    */
    void Logger::InitializeFilename(string configuration_file,
                                    string section_name)
    {
        VerdandiOps configuration(configuration_file);
        configuration.SetPrefix(section_name);
        configuration.Set("file", "", VERDANDI_LOG_FILENAME, file_name_);
        file_name_ = find_replace(file_name_, "%{D}", GenerateDate());
#if defined(VERDANDI_WITH_MPI)
        file_name_ = find_replace(file_name_, "%{rank}",
                                  to_str(MPI::COMM_WORLD.Get_rank()));
#else
        file_name_ = find_replace(file_name_, "%{rank}", "");
#endif
    }


    //! Initializes the name of the the log file.
    void Logger::InitializeLevel()
    {
        logging_level_ = default_logging_level;
    }


    //!  Initializes the name of the the log file.
    /*!
      \param[in] configuration_file path to the configuration file.
      \param[in] section_name the section in \a configuration_file where the
      configuration is to be read.
    */
    void Logger::InitializeLevel(string configuration_file,
                                 string section_name)
    {
        VerdandiOps configuration(configuration_file);
        configuration.SetPrefix(section_name);
        int logging_level = default_logging_level;
        configuration.Set("Logging_level", "", logging_level, logging_level_);
    }


    //! Initializes the command map.
    void Logger::InitializeCommand()
    {
        command_.insert(CommandMap::value_type("hline",
                                               &Logger::HlineCommand));
    }


    //! Empties the log file.
    void Logger::EmptyFile()
    {
        ofstream file(Logger::file_name_.c_str());
        if (file)
            file.close();
        else
            throw ErrorIO("WriteMessage",
                          "Cannot open file " + file_name_ + " ." );
    }


    //! Checks if the initialization was done or not;
    /*!
      \return The boolean is_active_.
    */
    bool Logger::CheckStatus()
    {
        if (!is_initialized_ && is_active_)
            Initialize();
        return is_active_;
    }


    //! Checks if the initialization was done or not;
    /*!
      \param[out] options the new default options.
      \return The boolean is_active_.
    */
    bool Logger::CheckStatus(int & options)
    {
        if (!is_initialized_ && is_active_)
        {
            Initialize();
            options = options_;
        }
        return is_active_;
    }


    //! Writes a message in the log file.
    /*!
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <class T>
    void Logger::WriteMessage(const T& object, string message, int options)
    {
        WriteMessage(object.GetName(), message, options);
    }


    //! Writes a message in the log file.
    /*!
      \param[in] object_name the name of the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <>
    void Logger::WriteMessage<string>(const string& object_name,
                                      string message, int options)
    {
        string object_name_parameter = object_name;

        if (options & uppercase_)
            object_name_parameter = upper_case(object_name_parameter);

        string result = FormatMessage(object_name_parameter, message);

        WriteMessage(result, options);
    }


    //! Writes a message in the log file.
    /*!
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    void Logger::WriteMessage(const char* object, string message, int options)
    {
        WriteMessage(string(object), message, options);
    }


    //! Writes a message in the log file.
    /*!
      \param[in] message the message to be written.
      \param[in] options options.
    */
    void Logger::WriteMessage(string message, int options)
    {

        if (options & stdout_)
            cout << message;

        if (options & file_)
        {
            ofstream file(Logger::file_name_.c_str(), fstream::app);
            if (file)
            {
                file << message;
                file.close();
            }
            else
                throw ErrorIO("WriteMessage",
                              "Cannot open file " + file_name_ + " ." );
        }
    }


    //! Format the string which will be written in the log file.
    /*!
      \param[in] object_name object name.
      \param[in] message message to be written.
      \return The formatted string.
    */
    string Logger::FormatMessage(string object_name, string message)
    {
        string result = object_name + ": ";

        // String of blank spaces below the object name (in case the message
        // has several lines).
        string header = "";
        for (unsigned int i = 0; i < object_name.size() + 2; i++)
            header += " ";

        unsigned int width = width_ - header.size(),
            read = 0, message_size = message.size();
        string tmp;
        while (read < message_size)
        {
            // The end of the message is small enough to be directly written.
            if (read + width > message_size)
            {
                tmp = message.substr(read);
                result += tmp + "\n";
                break;
            }

            tmp = message.substr(read, width);
            read += width;
            unsigned int tmp_size = tmp.size();
            // If the message was cut in the middle of a word, we return to
            // the last blank.
            while (tmp_size != 0 && tmp[tmp_size-1] != ' ')
            {
                tmp_size--;
                read--;
                tmp.resize(tmp_size);
            }
            // In this case, the size of the word is greater than 'width', so
            // we have to cut the word.
            if (tmp_size == 0)
            {
                tmp = message.substr(read, width);
                read += width;
            }

            result += tmp + "\n" + header;
        }

        return result;
    }


    //! Returns the current date.
    /*!
      \return The current date in format YYYY-MM-DD_hh-mm-ss.
    */
    string Logger::GenerateDate()
    {
        const int MAXLEN = 80;
        char s[MAXLEN];
        memset(s, '\0', MAXLEN);
        time_t t = time(0);
        strftime(s, MAXLEN, "%Y-%m-%d_%H-%M-%S", localtime(&t));
        return string(s);
    }


    /////////////////////
    // COMMAND METHODS //
    /////////////////////


    //! Executes the "hline" command.
    /*! This commands draws a full line with character \a parameter.
      \param[in] parameter character to be repeated.
      \param[in] options options.
    */
    void Logger::HlineCommand(string parameter, int options)
    {
        if (parameter == "" )
            parameter = "-";

        string result = "";
        while (result.size() < width_)
            result += parameter;
        result =  result.substr(0, width_) + "\n";

        WriteMessage(result, options);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_LOGGER_CXX
#endif

