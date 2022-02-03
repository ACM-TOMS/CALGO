
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00833_laguerre_vary3_mpx : public TestCaseReal
    {
    public:
      TestCase_case_00833_laguerre_vary3_mpx() { }
      virtual ~TestCase_case_00833_laguerre_vary3_mpx() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00833_laguerre_vary3_mpx");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));
          const e_float bk = +ef::quarter() + static_cast<INT32>(0   + (10 * k));
          data[k] = ef::laguerre(ak, bk, ef::euler_gamma() + (100 * k));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("3.4166346322604973042606006574787249780258419997787086036048027220995378977569555928351275027168301232740639201759725256624817525127027165762051686145105539021170618649570116301529987392137854084352753271702851711481321432574380019330919657193521414205623895059259341864440775608092645482503473390108522460667841593644641580216984427670775609157095834965345117503597002064325198821700985217802053454651e6"),
           e_float("1.8156727463532922283816790780652080072839748288051118825205516829524125145492225353226964168034812485328925949669869686176244485077204013352189188932303003278702292360555955509521827083588000599015239828200482160221055460053065369072500307459844886008027663574453249034740434180660284467529674700290793300224056904601891783345019345799636966512517702939148416646071487972418364387836318744419217641003e104"),
           e_float("1.6754462191319435841896004986174266509865259116219537896575268880964472265771170610068670568894402161268459199075754263692113010187315428246763942860682838014333624824443184277438381313436141663725560559662188871245282629097543209939426227608681883575300873503040615716205145165232493969285541245638711953114058757320719270216849971965171331986898284695683306916586872183626680669199288934761370887297e151"),
           e_float("6.8303931596600001362572501724922831501738814171936826182896438307741697437064050460325926162825402009593722711493306229709047143322502530403949734311777280651698623960139683452942692587732672353221775754714894844149219960835361692381696520666101093394828259131120967029889145199173585082773126461070449412804257065000037155576563900745975694140515482867253015911889337303527642814569188054983095764971e184"),
           e_float("1.2370953642636342023325958145161076054649639902104561090964213521173977240617210067122592393657544665658729940166112178735647130028195263024503410210524297847310125859835996772068547607299075518897555484881074397027794881975325334049948230442160146616369185368071695991663478276983916454433783625517301836761078206528993148103340259379901413309239505419454891871538021821372567872318414394870475630093e208"),
           e_float("2.2376072947740851329002555478427562935826497399124741932326305370547804984357922646834317173539734747502073725242149112720076535924791737883539311137740491093708436480797395923553488266159547083289987646134834620843820863867640984598250812657992946023758257666815240378647871160178392299504724105788620345245679412149057985245564989760205326396609709530700010505672580436765359606007524199121020461458e214"),
           e_float("6.9214528488353522580897263487393595092209660138272929439218063742254642687120710899394151426978233808329464041276297139013650844925781526251325097857355272091131926264821562234207545307097215168566726943383984500283901612762004977599955775835912048960534829240111504881143252414098524317357619541244354029010301523810617640055771114392992709028077159490545553105697567195520775014460151598810298122877e219"),
           e_float("1.9869359369005751143438637664269184828142366134444964154151988464384949447284528020791850797331102260081638697086427595029528330417300400835148638956912094640105993590061886040414097611164311311328227465955535999806706691578780542016906414508028407775202669587375761511664529995648769646144271802641616489373194575131081519239142617834714845068863449232001143590871541657337382619170344771413692309358e234"),
           e_float("2.3803536278221625029310600207412014076232260475048602093843253137117307706002443976317796263102941348235405575487342751008972811512593222499390073909638135446647407024340715720747512559139572947235682946063324413186026346981384507866147499426827193809581454479985639859699462908823333828508158245559125411171003728951562570691520277359868456883822093313767175952988077518669994602062784545805457789218e251"),
           e_float("6.3452329347307701821105175445547039633568970542888008430268761036773243132288993473323806358859865986247736641251296612466606951545492377995547088667516347068328498149649119300347666593504330585570530725563981458875686977392597322621583468252733322623726748220725327279701298777544642889469862178573872012855942020975806065689881646139282857141207712458835620343249065744420305375595925789844267897341e269"),
           e_float("1.1072262431125404602766676006012620806281816156997872108545803683081355505195873149812359722686653927075990524283750881700255910013155373354651521632963233507564107555043196260136603100937521533685528598279615891309895617107052092279881246540035101045450181925393647709832476404551533298548739312303675238432881054439739115513028385669488854167120795470604709740344211901551208837432218287342231120834e289"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00833_laguerre_vary3_mpx(const bool b_write_output)
    {
      return TestCase_case_00833_laguerre_vary3_mpx().execute(b_write_output);
    }
  }
}