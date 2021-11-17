
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00601_legendre_pvu_vary_01 : public TestCaseReal
    {
    public:
      TestCase_case_00601_legendre_pvu_vary_01() { }
      virtual ~TestCase_case_00601_legendre_pvu_vary_01() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00601_legendre_pvu_vary_01");
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
          const e_float v  = +((10 - k) * 13) + sqrt_1_3;
          const e_float u  = +((10 - k) * 11) + sqrt_1_5;
          data[k] = ef::legendre_p(v, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-8.7230876733645892179509818358063268752208266881890689024837436348144181227422386769479337294593225106766967285139076268379153295031309392267030046259551356315361655097728663603828456051636078703195809403394132039627744957208665569080498976428592092709497598386846494733874556568105374379785971161523195208277786438352316077947721055713764795134835671333575174580679182679552586051235403826322405479444e316"),
           e_float("9.7927525506474878834655327303914533402644954120817572613469431403066811842812959018129171473806505501392641936398262834586187865954656123974077590607995572976748222732311214953285490018389550607838569956211996310166838795958408356924394294146697028997762756982948028926207013036340341525322839762841212352422368245063734225729268646608208036343756215983834445272039856850202749014380355817750491062049e207"),
           e_float("-3.3610073276470566891526730131061979903152554426316318839559469309752139579412437918312970696149743663092071148128948997501735656660178753294811057013704759283319773121474037419524322353843064882126401746427261113637702181394701506680846848259628424579676374228427352368936765923836137364857919959056433208726658402834285851787175222738906126168280017939726813195581826184311355276099228447865093106446e172"),
           e_float("1.3419019623483117665149029785597083546413891986252991361249904914306149468264707805443758900778873564498205174287462073385542210753902796429709456034306187364100270313155733470989870990911516274483221462388638008717542813237885677117786398251424128740017078404831831986302655001872837714923050456043726661492244003029931996937482145272617050489354146917484836555969139924630349296174027738647119885212e145"),
           e_float("-4.2224972375677729301630059128183375182387732331255454697093510020301319578163791809777856537056998972943021202042591708257474821874893316784673736690080945880326069140913072575257776411271304195795006623020651504675758750811993581403687530263875496549302801454403409264503426364894813044263691687700337978901465120672012756501618402324130206455875058169556825595551804048516163898610922357298867316074e120"),
           e_float("3.0409213168255910130648765361421415067827475751226864318473952993323229074147368013651250912803747447377262773043083953698540962147724127852993175021005357456254681007957527659704135942627022771387573076285641908678378555820812536352351197699115480215065671964769407697440827844650309132112091889404357851295764674295939134407614879138034073296553967195664699502063862402874720901165273631043121300191e95"),
           e_float("-4.6015327535519105235526558309774901368058274166250028061524939138387689408302055442679680844382850814528692526384618421612977863732667164759627053972324462787601583938613268413026987556352730512427431079900663476064198354759090871176780910448269581177056922649520790740112580049799060287389211912530647651592406576789149778680299268835737694315541024089262353493074667126818977389064202845744591559937e72"),
           e_float("4.2509530655060990058894578893259172134514072588336334831907985632486351271884193753385729373217237115451807466292674575282697872172714994735492351872837925757028506229590617273122305870413332301638370751576094932080274227174668712446891973455812183712243267669144812474496061004137836126602118933385277926175059627352829896589511887849949912686851388393479489172053702877021338662270079518877647621192e50"),
           e_float("1.9945880938946018294276538056839254052603731645730796745200349572771305592646411277719500343952689077199775351248486002233122752159510880849725491681026248209560243377156517146851387660981397777099395491732043577864477845536186010769913110025131594187427861094266733332055496084828906754159398320711028081943466359729708249454777422832725259527852021740421868193358909704899340918032181057685769646253e30"),
           e_float("-7.10913560229582249236279165226904821810256628701376480973831057270137368084459629743677591385013930517908366504213613538508071958239104610249045075953711816328888569694836610482774902288002885715171219698490891416133358530484417241019566667720388142196848915589051232852985953822023696491460555222424138383156705382658787666737591202722663505782980480220338949043886685519131085739380059284312553278e12"),
           e_float("2.2864339137956973572965932101831658897128659236302833052295452721525158227656691911865327925213601170478121108155801304010852977265617924867919809944976478876513887519643117891067262638129348936681063606672000421496634690162129382845936244098782016705861431055760994999464085885936042127273893633600681893108981140878536496140247830798133082904723263597478465392677377617121450535738147438877085480356"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00601_legendre_pvu_vary_01(const bool b_write_output)
    {
      return TestCase_case_00601_legendre_pvu_vary_01().execute(b_write_output);
    }
  }
}
