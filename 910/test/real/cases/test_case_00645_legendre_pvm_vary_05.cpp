
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00645_legendre_pvm_vary_05 : public TestCaseReal
    {
    public:
      TestCase_case_00645_legendre_pvm_vary_05() { }
      virtual ~TestCase_case_00645_legendre_pvm_vary_05() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00645_legendre_pvm_vary_05");
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
          const e_float v  = +127 + sqrt_1_3;
          const INT32   m  = +127;
          data[k] = ef::legendre_p(v, m, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("1.6565072344110981214807224507395291876064276711622469834723112233364431353102824332235556854584206541517606753694330540802698039315513590579917719823226538700636227764129775539978189416541721700274752354073053240019371689108378006615257473902916012914904866263431658549808011920709867789823644479807182499756604837209139850686908833925682349313890767668557500406284719955016482000246465671577893082656e372"),
           e_float("5.4258192252944467720022851845149974031309644314743527797926258401240979618616873462779858387341133602978960184470659986154316579447100629568681531971871363436272401607902930881200501643321541862024137337646281753887485300792257969037710119299958536263945554560552088488066305453840143156909052229869889001395843510287015250168126957973515609291304371948530192338896340367618786223771731990623193876055e276"),
           e_float("2.9587467766680468822719731978248325891670791683212310452778367855423886544137736460446741313221082855730563226432164020871876001587806846969880991308419316093609588009702376800210307860749259037686094018995330775628778385424004477877950098848105978060964904554086239645244691291356614926241479520284725996909834983458868307054590806915105407522928596994707629229591982104961520839726671108285360065494e261"),
           e_float("2.5833813503403828540219647349538243168358131417134666800651377908336385786879238347013310282247707904073481943095821868431108444743838196335989841930344030691594792857640905431904003961792468895358421696787292527185015715852388225039157062154761817056539741836159369276814064230834103303511506693912586718004960612163531556055773882960760398784755290185379034056922993230044863951863065012800177132783e254"),
           e_float("2.6523501472455506898569941559421738634442653385904165488536688968634428422696404678761597512211242710806978742556740349524434967550497809784590877086105373216444504836111508162010227247644580681331351200342512281881223659127816381828310588060393745765150879134840084818511517870933425982549578170401051850498122346056568830439444859959282457351647658816792206057291731393549022247410486906883872031683e251"),
           e_float("-7.9862756438886187152242472958628740854246949923332233409352389242771782075051310507765350086775467899588550755187834828419847506295406826799795917843914064976162836897884685521026567918164763691812435733332311509186201253696609635021086353392354046883340496531154834034694652005506067407160606405578376790144516116884476127467224671052179661840637812247225788657960577383805752959860007278264487338922e250"),
           e_float("-2.4077446046054554375046298736004252138952328214758933831639946418024888258723235510569424089845121254896648977609492306491190901370652659618282379325139066849360597999396172929977429923358147284823125299767219545785186044221572547641754665435275862027099209638342316590881893692584887800068171701516579079919097101549862240572563977770672530118776983998971213091647233672696275600271512801511404257906e250"),
           e_float("-9.0708641435487252953613659117532821945975524144407411669215972810727805006730340815649605854407901916781924067909441890058867994819909453461534718436067465598050666596168636831620869014031281718417704751115823109453333571220231932751512317211856962325622591408442277988670425383176170600953398219253727347779514760480857014682319553705404744163872408119131441405089965050083292591799150201494870527067e246"),
           e_float("-5.079798257662248950815986847797677298670364321469260721186756819698925301625497880312301839440440540445912353644093411059388260656051061827291518906117531422355960685224291522473698927654033734460163653053699015171985409955293256813652203658267526831573671472856671940414418570792257701799711376029719348887863866376569114510940708788717888598240063831516376167074412115582306708872960475707779068974e239"),
           e_float("-2.0499193059477495470950171974251296983811206933608517330167466999783667322327383184999346356568447991995953788579270477723659119140143003515932852000391718540875036509973224057078755875171260913190526796901087452509891873855344244005165672788343891221880230975431259005382903098755746892413887565736475416015908929489178504909014414950413856180327708579619829168117733327145608717380090204771593630061e224"),
           e_float("-5.3373400691512077695287866984139868221390559124593750580204322407331683383719151333213239074714091511577703579086012856315424735994339161820208039132563425684711016780689719039587378112118819884159607047132981746155732851133095285602321945308322375700622788444148336316665943873722078814356638337153397028988009876667984545009327779677785682306971992841316409587098065152249547001910393423955200336954e128"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00645_legendre_pvm_vary_05(const bool b_write_output)
    {
      return TestCase_case_00645_legendre_pvm_vary_05().execute(b_write_output);
    }
  }
}
