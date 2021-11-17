
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00611_legendre_pvu_vary_11 : public TestCaseReal
    {
    public:
      TestCase_case_00611_legendre_pvu_vary_11() { }
      virtual ~TestCase_case_00611_legendre_pvu_vary_11() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00611_legendre_pvu_vary_11");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        static const e_float sqrt_1_3 = ef::sqrt(ef::third());
        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const e_float v  = +127 + sqrt_1_3;
          const e_float u  = -17  - sqrt_1_5;
          data[k] = ef::legendre_p(v, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-1.0838318554954524895755795254387636499808485881433717427999015432044336788174655720456617777388903462900946580576125827080276294753565132518149463946421966790027105234014788448623531099397671549692686950844074385543033061111999146162006423704881201766851095546598763983439882000199651822853671764580413658942918363866182461207150863517051210511612375900326253209189631646371914243185501052378968205469e-37"),
           e_float("1.6417564845720513127616729414463446987088339076298615423923558392766438816915874791923611408357298906049849731141936837979608434352367887438827881616980703634291620167513593879996305122638952996551657630053545478432813835151660074289605497444297839332101213442343827485661199169003211562780522191005441369650961760300816145231069384142049124440543646396730504607124796414811819333066487785290794517094e-38"),
           e_float("-9.2430096456504316154897506491115028202601804550670086271345847235206281358170245721900138680189126045189488302643588932834617627719064021034355058501435340755858007323312416397131673324757653630206743871438723636242670094816368943085519030471664294527089983206978610858584745504356579124962392464559899092633729827983526706835008005006289838690076490334967050894887434420625880221860318103568936373856e-39"),
           e_float("-3.1971392679106007282708192065744313630584119560351772248140135581438210769128782545738797974463084775182629010477194069671990670077800026129748821835444981643179396363990956202822726395903042563098981680878270835108079009629883696214831995349966837302896467413179420736629730774697834175232702819536865103182980367650400356133405771182583828563795971136848556459542079421206223060207147645989991179276e-39"),
           e_float("-1.2731587327687953860923952424627190522172689172373211627261575941921290626037848705080479508881256065858107459780162371760160678134074160583709283704318847780374861955047867385733725450260907055833935817480975162325939270007580665971128833593717628301655980227988968081134559885184026104848041069490022103846896508198924995034530100425567844294837308469864506391557023553683011319872291098988385992245e-38"),
           e_float("-1.0982281011957612444979067744809548425369353323266547750867376032628296966092608676964107357958400635515175817610284532867025638848755064721533955653983735051651453007670324765608116486995232186239047628847221841527321775159086571810957041300401223748125720742445127767578656092775856364459342230585309137568968877653503932119881770166151234832594941463823350244067141514285110640937766817149182535632e-38"),
           e_float("-1.1024607976574542167196839142085624652348857018472935313210692402228531625950831521240820447509272550466489253331903997590298627559203205212591571967689652068885298632833183680232300282181786841110486172251817688536513847533270033580492403678635980814488087144335735211744537840355382111350241855484886993989913153337102634583066894884177634774700128465414597588241999346500051196346063510649409705244e-38"),
           e_float("-8.0549948823568384406506106882462475359112793890346190086470373676526466871238437677197125039308734517327360342834181057847762642157314128282668468873944849087372038795829245051197460633222325209834197914216748734125239830336353589777935885743247190944570038138899432284749062248323016655023577089683147303576287315744665204219102552962728099327060287790424499047589416002402626829707240136419413324975e-39"),
           e_float("-4.1890574303989776121198551985925414265289291674719392442796504864660374221971950334932405442439435592785082216266477789212619688168993613286465730268887700448337862038264771221267480461913367918265171952156445034122504508531321511284677063876273084494290831374019017979579730493204766954735221621094323337092602796330675705915467976862328573759532846537606041452101192681326921596535957483554025619203e-39"),
           e_float("1.4693775400934274407844507746744940855772449889398824971880589264734870478332596550511134644564427342528203810948930372588044871371484775257344063059516941020867189771525892632067163641434496792028790347679462194354554739480080029348642685848856809804356440954189223592790979045618160252022171783369995529188530070097592893776646238959037506526166704606885808513463277218132646801804518113381390343839e-38"),
           e_float("3.483112700809213386893221803093352022442398816380281787565417817337859455382404912702552950721760142601975572651850901236151806961995737859085828897380036165979743872678786071224982441003111563545809352113843622451730300389644705954863305758289089234721014960128076770291163890955962709538654918370490066559466058747994398368857376556742327243307328776822627897043490309910277146402423464867312085295e-39"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00611_legendre_pvu_vary_11(const bool b_write_output)
    {
      return TestCase_case_00611_legendre_pvu_vary_11().execute(b_write_output);
    }
  }
}
