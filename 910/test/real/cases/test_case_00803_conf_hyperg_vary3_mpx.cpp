
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00803_conf_hyperg_vary3_mpx : public TestCaseReal
    {
    public:
      TestCase_case_00803_conf_hyperg_vary3_mpx() { }
      virtual ~TestCase_case_00803_conf_hyperg_vary3_mpx() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00803_conf_hyperg_vary3_mpx");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float ak = -ef::third()   - static_cast<INT32>(100 - (10 * k));
          const e_float bk = +ef::quarter() + static_cast<INT32>(0   + (10 * k));
          data[k] = ef::conf_hyperg(ak, bk, ef::euler_gamma() + (100 * k));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-4.5238670687161214455465152113574345123566767001371854118437610085914242583289228017225544048694169757174297959240873072051782947159434991582358200191950181246213789231293445356072698198509587446420865185638718130578682344644409936610868971186982877404590873657646749039238793910795357219129416586649563537688269963978973872986924715149745855159662706536849524458923922317603260647606981624365636524515"),
           e_float("-6.2538617452797367671830972565568381364069408439102328271094617202105240433942183354157208305074640690036596695512752078454051706102579226072230809809904817803340621732115371624982327398456241206816211442052965808119669707719550028986629421894652665940453537644092345610427364825376043507245087191906765990963580086411735859858517234177913328554233033985024253470979135958127337429027644157618085368756e7"),
           e_float("-3.9155835644934078328526944649777487854905937324964564503315141561228769414589587082384662100113427399633399209319473385504962523924543787295923053765283159965490434557455669654930705160734283346858606334135585415450088730343843439057905428554944948099791556001304706623138936083156457162843056761754202017627353824199608260653482565273432255106996222272642750438156105988966175431479236487939107175316e18"),
           e_float("-9.0263503068540164497934384306405914265213187146721424954804073381415560841811468923997413491317422526226641253926527238683512520477570508167357668330229016485663709756088765069909560011728572332889490961530721514234355696447125033598955453323820203598901099729435234555925970024644273873107119092434463484485733425978456666023638497425638330712239536351840255028910522411225708279289710106358241995988e30"),
           e_float("-5.7728006678407635349803538471842001606437737468698117088984347155669230861942509054889289986123141195163626679563514682652113494595482362807134573115414105271369793451625152489563133177626573908268162778925382595514617234935742761436520139852538276553702287530947201492043821931052070557624428287424903596856639086712879272177429151407978147301518454863137683401917906690187856881722129958471641903475e49"),
           e_float("-9.9370820048176679778284968847230508615429040428933353907494123502168980832173073870164633271307687831598381054453791264484552227498397411157353056141759531993259441399213923816813279310194660569560784873369926635985177634541965634702953471105674654895914754000920680457201476656440031009719844169735198237703973416814007805523381744128566152519271964777033551622218913045647434131337738299451925517711e78"),
           e_float("-1.9050314941404699262573850460059664831418496545923952283037279124917178839247365839266496936870849345132714179748504036366439475093770919411153087120566518853938607783718807947680881698004701290698355842877155967963096255313403450847555210753925232239498363728588911686554216615650112261742648318438436129270281911556652207264725193414474043632432974523620074728817909515989094126215544097148889740996e113"),
           e_float("-1.8988121479219410611950768998335538365833039061145397132714905440047517235780832353337475553880660990536890856485625024109090405727129494705163114432443073904400301462066596170379680802215316179831369394288021474942559223588953985333608552059580224494223781530744261198707576813931478760311916614267125421475344974220532312633627562405649028587585226690482269978838271597542347127681797901671017177249e151"),
           e_float("-4.0727663578323767828695769285179901010573374207922541500317053397016286618918390790153935631049993198566658178425300842561448071784252806632358900725822744387860316210218454479586173352228257777999033475270480656797701379801028316643188634008791541656642773160504531056203135877753717561758058358873687592013910559455068628971867484804627746489058214245058126436156070233487756232770336383648522603458e192"),
           e_float("-3.7720146008975071803458592724371803022577159536765588240552926436351694540583327589044132249236954038442550023785143748156869615387350456051557179819023812745903231713677315262474440477946783400940855954740844920330869224357640351842340876973788764162736387207764808982076799455950283955003367148841182043852375347484869584691536374342817477697755629315917027347946760763927883864921904143004243039447e237"),
           e_float("-4.9227409711113531906810259680833882596202535330079853503784069640318335597190992239002321497600693473270903706997476731861454782578802449535216104967234690533977064323263462396088277982050734347967313534099674710295537363577729292929823889039299829478246147466401710751246197475494981597479161316627683397691478833605730601635293161157627994212219468083832738470779980396614897947769524222145807468123e288"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00803_conf_hyperg_vary3_mpx(const bool b_write_output)
    {
      return TestCase_case_00803_conf_hyperg_vary3_mpx().execute(b_write_output);
    }
  }
}
