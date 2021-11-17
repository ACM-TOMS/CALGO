
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00826_various_conf_hyperg : public TestCaseReal
    {
    public:
      TestCase_case_00826_various_conf_hyperg() { }
      virtual ~TestCase_case_00826_various_conf_hyperg() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00826_various_conf_hyperg");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.clear();
        data.push_back(ef::conf_hyperg(ef::zero(), ef::pi(), ef::third()));
        data.push_back(ef::conf_hyperg(ef::one(), ef::pi(), ef::third()));
        data.push_back(ef::conf_hyperg(ef::twenty(), ef::pi(), ef::third()));
        data.push_back(ef::conf_hyperg(ef::three(), ef::pi() + (ef::five() / 4), ef::thirty() + ef::two_third()));
        data.push_back(ef::conf_hyperg(ef::twenty(), ef::thirty(), ef::ten() + ef::pi()));
        data.push_back(ef::conf_hyperg(ef::ten() + ef::third(), ef::ten() + ef::third(), ef::hundred() + ef::pi()));
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 6u> a =
        {{
           e_float("1."),
           e_float("1.1152280911295835935556221710652859408561228560835225917867044262566651717588326053116653022313301769575567737738229808189671351357799759160905948859113471574935893526705551753593931616020063314966324388328743758767384913500354918510544317732785141809375943319757539270693272641042365828367838323315208073715337605134182267409366655530889751838572253929394390410556331709089721772803761044414476391441"),
           e_float("6.1053904278276446991908894668756736986882819680831771639906998084411542128907562012504047426786539166694374741114389367417651483759512445332877788304976443999719930476989440957063698383850944623237095952504972705833099502514609193647061784640139494654928858655403423908052752745656968247196691712720950081993022760465066659344094051696056495923643322388981398169589369399884963234219146045079745315876"),
           e_float("8.1238708294703393503755429136736895407956397104466463581971919898835456841455313116386303739772918691737633389528277383377739886169016039399246275547321476578907083811451374758863765021589005914026262088838241970018212632090799243612382948756662887253885194346439246926593387746999845920742090554760144524195382254230882528797783729379961216395052236306525833502235969786257990435112189894394611729206e11"),
           e_float("11163.141080269263241015054107202083061306962581112136089781111823873605171228729621892725393528803299459509198813942646047349493930050866144956594992826001684008447772312805296921140742715673553993892713158499270126838434616875462172440046847520341678583493117413888390330624438346394686979843169323990147058778036538595635695066944946723491671811189487132514261271391753400012976650436811825030716758"),
           e_float("6.2204892539672311042182433524599248371305688864117709488184437790602799086817528027435066102572536604508115176291687726839103499656104252381870000959864283768239297502712248229633292225879641078999435252740167979522098151549006048050352568659150489201252308170088706597897666660813246613267604927043355019563020517274461518602405169879143600425690883093000875141709405130053958728259785808134544906554e44"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00826_various_conf_hyperg(const bool b_write_output)
    {
      return TestCase_case_00826_various_conf_hyperg().execute(b_write_output);
    }
  }
}
