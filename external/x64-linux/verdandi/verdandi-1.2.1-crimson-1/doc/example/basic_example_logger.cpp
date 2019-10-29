#define VERDANDI_LOG_IS_ACTIVE false
#define VERDANDI_LOG_FILENAME "verdandi-%{D}.log"

#include "Verdandi.hxx"

using namespace Verdandi;

#include "Logger.cxx"


//! This class is a test class to run an example.
class ClassTest
{
public:
    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string GetName() const
    {
        return "ClassTest";
    }


    //! Calls the logger with an "ok" message.
    void MemberFunction()
    {
        Logger::Log(*this, "ok");
    }
};


int main(int argc, char** argv)
{
    TRY;

    Logger::Log<5>("TEST 1", "ok");

    Logger::Activate();

    Logger::Log<5>("TEST 2", "ok");

    Logger::SetOption(Logger::stdout_ | Logger::file_, true);

    Logger::Log<5>("TEST 3", "ok");

    Logger::Log<-5>("TEST 4", "ok");

    Logger::Command("hline", "-", Logger::file_);

    Logger::InitializeOptions();

    ClassTest test;
    test.MemberFunction();

    END;

    return 0;
}
