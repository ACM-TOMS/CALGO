
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00613_legendre_pvu_vary_13 : public TestCaseReal
    {
    public:
      TestCase_case_00613_legendre_pvu_vary_13() { }
      virtual ~TestCase_case_00613_legendre_pvu_vary_13() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00613_legendre_pvu_vary_13");
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
          const e_float v  = +17  + sqrt_1_3;
          const e_float u  = +127 + sqrt_1_5;
          data[k] = ef::legendre_p(v, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("3.6901959631345461357160197976522575764753005516901970090999415712248582585921452068789841611434589337413144976251674418066006778371178262201871135599085543551911210041426537743291748744467328318404956032730118270044706609690533017159169314049887213849635458106340360265354330039180317454891883728915222026169683898849510648078339518078456994415325401010853785535835461583834353849602992319503310351548e373"),
           e_float("7.1180235286312385677288474069852843287889176214340235910511666489826607906454792687700922684892044842043166911164254856987409118439740315007355922636551203235466076912816164016252193766005654600646188882487632765975835417051594553858246159280910582996070389881353531584841194247758911106980458170295869241930079621646147629058356684774077887655577391643260202332806035288715556249442885045852404947035e271"),
           e_float("7.9793641199662906943014299707797375458243744594910977432470579591779436077858986836526445073304576238637645823501027778305086665129972379200093579595379844448926217275227839568660661697624914966286865795796599127062827659771013062064042577493459928313676463517247102493742680919944449616292679246401114672978977557506819579654499404668714760038001213564628439261971208721240069835071426957865675506235e249"),
           e_float("1.6394433371826182460059207557056302136144533979334182572384372937204402363424755127796433843851329718204395070425715101363336576495585060360099484311822372378349203052454933500185547313602336463699924847974143662667564407956915984793618119448793685784646380454944918103654104063410115966829759963944376269578685325989046327074950297810619921661045835282323952285932998239243732394287177961891155256115e235"),
           e_float("1.3996080843774822931009998779332766375646793950193678482384145887673598781403242507264377573188784532040568557408152890198717205595612334045664507137245475998863189043494970379943444777909777687970504814363655777365594628032333623892231423766502512757527597749081476840320625955230652076231899834803623876386476881027973097445400606535656213516305770845806131827805160991073226721862258573751565164505e223"),
           e_float("-3.766929434793202752550884801449850825543682866074357241962070330975169710728878400514080761398415987022116830495223646833399680988485123658807343202531975122520204192073295782888255757662554858527499354117688146628112973253016485977766405886106581950143192787379018657004672851211443000569686901102418458816690044757822798065897084640671429710774948671349115003341987173824678476853347571189270687957e212"),
           e_float("-1.4221913549420363945807257663142029685533182349469659371305152329754969491392545751154127516946339106696633499548422859653767334571601006290135084218264574379578652995790656135915172848603628111043891320092611614733276085366886136755000845622374425768327603181832124304361197773895587410326654957898725595677535345747186887387633000586936116077004056931599249862361854797046790975124282570002753582803e223"),
           e_float("-1.6658964513594473254722707837123339160761472097034904275880630552363330440598842577813020668719967000399654318362420217474971307010975611532787074397517459089361723812126575560118038267088657779243611913659620194827439299714641604003711731155568409153897606593163352073583458606667492899528655309122012911064673877076228359193107608382859760903343242818785689525014478666771610056913555435865220331105e235"),
           e_float("-8.1081145472221058788160908703400624253297270109661253418105487211727476521385379668354948886681818233276060048652882342182951873083351734620143545351973954954918377681731165894305311341607524271159435074429730399410000468270705171994693956301840438885077733593301884233720559621808569832569824513424312946048547107925380568283123113193568107365032704482930581300343041565036559181709899003535009899844e249"),
           e_float("-7.2328758598132486933467849005738005401319183985289752214893029544764517452834277350153995087068544895104499861821465966746260883020057362075699628512678830059571185836746095348804783105314362975495967430275033886120785329747529414234372799917055162427658903082109624658898461903302980873201975527674967571114947873475954414041543787733001515480027467004865656439904198720128530433993224400563139369853e271"),
           e_float("-3.7497388414602020381962302764400255627506637075897686906839619939074265551870334229180300585561213732325948305647950731674458319034864104178157750867239211433582015019593342230838504202925604503691600357212280198625209541708024229180088436108984486636575924753695646228109914870033997828006286842512730381418034049438862821275804665396223157465972579928405854408259100900961883396945807581231825176694e373"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00613_legendre_pvu_vary_13(const bool b_write_output)
    {
      return TestCase_case_00613_legendre_pvu_vary_13().execute(b_write_output);
    }
  }
}
