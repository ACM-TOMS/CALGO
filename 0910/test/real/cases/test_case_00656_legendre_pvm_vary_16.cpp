
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00656_legendre_pvm_vary_16 : public TestCaseReal
    {
    public:
      TestCase_case_00656_legendre_pvm_vary_16() { }
      virtual ~TestCase_case_00656_legendre_pvm_vary_16() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00656_legendre_pvm_vary_16");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        static const e_float sqrt_1_3 = ef::sqrt(ef::third());
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const e_float v  = -17  - sqrt_1_3;
          const INT32   m  = +127;
          data[k] = ef::legendre_p(v, m, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-1.1446497738474759697745878591939280177703684157898046544299353354360764486017220108309323662268785762654029650169103121946989379428987198169419218945188407539025621200904742340188954174911743950856981998915139345186342011598389492856901178660079609533438068692642818554690146458668361330647329039684093908206862140089747413852797182220023602338001649530079586257407403321613421940664635153378809678495e372"),
           e_float("-4.892856997023379550814743925414704161133621442034882774663643959389085919047718116028191296524292131388521192888538872945926383258354911724209006956486744723739691423587918284939372417323738508765937303503144014060601184000008680168724587680863847392020851483258552661640694066513988902965061526562984758699614421745102650587560547759095385550149046050334760225857764701966924873003663568651612446545e270"),
           e_float("-6.3802824089774330441043624548314871473905665957387532153268859585234136265200182432449430059417968288448378189455948933659204852526887409294659293518290502563482367460152085188067010531830341906223383059503878813139192878777246821876637956594936149219279117941054826961201489299895633591652150073666199856379574639436333301798091566392577742567005141149032323718588195103517892204530935801226864573973e248"),
           e_float("-1.4379424153976762349271950176582276891291550608533930262859899881894938205671330986673591616667118677316079376828855086258333136199103554344003631357580346097309274226098466274608796629942686198517493078368854831392691756613056554070556438810636390807655734353551106875453195005640594967151134444966823059954243469547847679086387099928781275492247374965731801701641650258011960335511506153359936755756e234"),
           e_float("-1.3183679479119553992566680457814069365975539614029042406047281219391922821774370041843638147880242827803952148116844700527267325611344661640962450560577144780778719053857098759956015633509132507654961417674257395335901819611784575965064869992909106131374499820892163209729063475938920334957842371472043489314743603320805050421770831134504609426680881283608759686613527387889884983005369758003180776388e222"),
           e_float("-1.1253129958128804445532569480271552234031322935923931313360146285759263015708078563320727051459487055774064114993680744373047750843190776223197004452327831826650120373627882774215205258199414153985932738478266168659506163177121602565682255026716838026476990279871411536537780570475404778936243155291353969031060566329900226632385063327488225438890820057745659019726983826927071542542081449419069441373e211"),
           e_float("-4.1050664427052547056814869380310460624822193495472441764702111493954101165998517587756480349786926394071075133103759049059031559941342974310758158739500494722851716226921278256288738672484113999619186763411739987633433828582690139506545527796501736125483821598817321506686556786396043040920284421324539664232821005168556546922283864601561947268390574822202400613823822168678833411870193382883015166715e200"),
           e_float("-3.7596379410555002729861493669374079867657747922095101212205815605907294199670280951725197746847983715514547552440943869463568161099799414509586944164225726942640678872256101680862312858155666851011838224380903240662147981281823008242476400563837843006104658479241980723621735148551674240154351487194131445272185751583841431113817349557018994102046433536688412342637912441871466011764477226251438574718e188"),
           e_float("-8.4578805688200085334241322968437055257814123285188880060931315251123618887191427623935238673810464456500711417407287827801497414730606011625314698331648924063478214978027282452793613214948015020852434256551171506557096174145720896554627612181770143590951341831638523031526296979180484589766219702231480587160315050974832561993344130833794587816535035327387847977774521583280975743708009933196118696068e173"),
           e_float("-1.1001185670284871726635897031798751008802256461628963618278250833246013116154857166371875452957964841865633007527471977238844385603534983056574318823757173630172732042128395193310046474734385775646311639271538358464994358267838080427003326489116290022270780977380751479997338554233972773145675610435385247258818112284242456120541836550312810197099931989244632698715249309844109912020059091779152960128e152"),
           e_float("-4.6872842599887632058393067430978338568579127753610829646129109175176534421314688240656274652769004095770651408909095611432944770778808011410890965689995750616373933315996904721384142661689882622784876190861956308166286310331139510731796252078273236571030437279690622635114024865945196610466486539001830347065051841771715653497809383091368832125219921153313469325518237015520975680720668635131834803843e50"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00656_legendre_pvm_vary_16(const bool b_write_output)
    {
      return TestCase_case_00656_legendre_pvm_vary_16().execute(b_write_output);
    }
  }
}