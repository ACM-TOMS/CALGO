
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00512_bessel_jv_order1_v : public TestCaseReal
    {
    public:
      TestCase_case_00512_bessel_jv_order1_v() { }
      virtual ~TestCase_case_00512_bessel_jv_order1_v() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00512_bessel_jv_order1_v");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const INT32 k2 = k * k;
          const INT32 k3 = k * k2;
          const INT32 k7 = (k2 * k2) * k3;
          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 10, k7 + ef::pi());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("6.6830581431796544884424968404019563838731737594663538567396331103387527750677758240189561173294702653956400547852844837553803435112448160951427678652622886911914190103131417177665796518528819225914486166232518063553738264421475684951907810944043451356288487755544024437984294664541443003737747205104649106400925852665486828970395404991612050536473460128559073188103307701853791754756247119086394496628e-6"),
           e_float("0.00010579493217986162568278551075626351075628299033092444916935433972188187941827540347551693077007039009295871007676932012257829460375894562657728111273613310011347670541346094583651250929549970135868249785421976281351502250857394823259402051629345248208057858862472783864612061768615441708156074233789475795301736615827741878496583960092817257232198592506845866805572376640449548843594139088058610087632"),
           e_float("0.033505089868421750732399926476969306904872751765476354609851902940287924191320635012840329068761138807207279511402547962052568529751313700959528157734544293257767947267610916941982057555451940266234630764692058761936970940429249493810498721663914410730251129338645104364731485714738069100345276323170368911869065296648477767219327309320013961785028838198265698108943428105785518026025013151172510396092"),
           e_float("0.0059374746937073024217451044773360438258850397781830041035207495022307206054409645473086618895067560856647863363139484971106065721381036691611487692654130708500829696017210410291843528715959366785618840959991820601141366555726884563885909834824543599822742274437376845687188460639897465184296225226820453033463776193277154696493906441245962319354691360681419398701924004752979153974180716574085413224205"),
           e_float("-0.0028584428673578864525155789292345523828147946273585484303183969996293887166581692611292553029243456635627716461980957234537225501410111922943354699733750314297766309447410387418256436579269764414645677276857140798418704872074123851460624226018527484980651739356892502642928720339255757534988847375642805214643948549786471309108661237384707293267514110571252923566137572344080563735184490025844954600471"),
           e_float("-0.00069705276508763717203313132400016047534279977673099608294493573181427030996817667228594937211639202169786187472080960308300849510852630860991775073653185563816538490182649250440959515029986322119670860748702448531148167662796391492863869121955636936803947071950843913191599018229464158468623989241086242942160574792497869291914573917467167865420939773377338135452120016366295422674848205534690693973877"),
           e_float("0.0013599371724983568687654036359511077991705484607873090227576846320377821855278039847741860580153749466940268130441956757979663510648517158110024686047353905901599766943243440555975269772335817182410369936211651479079814546528709808310724200596890671195270797151645005941984906659872050996585787540700153490588116965598849301007483057985973609803149931367698969402531649641691127156229477938198383716966"),
           e_float("-0.00042353893188743343878158169759915611284362827741768416091834202079169167050911395198780912922727282980999055362743952546259080269478719373750004330772868501717902281501580086589686930559003069254829643668783895090891780609134727461085546667150973760335193512262828293959836613742541801889570984749319938944781453942114848149332616951053620538012078596960969264733734130061326888872816212590810184888011"),
           e_float("0.00028910400509060439500621071174234573807306570115921911602600906708372829498235761625288165085804614884011529177267035225950746309679278968401722693881228160205281907526876871439182132840500185401231342064846836221637328236179709958253809046638347545325820745341124133540663692563471051479870552017968368509463099351501786764576847696306854236441537406062265621736365538299759204071688570347647796881562"),
           e_float("0.00028067262685083573175604853887382727868563253766102153980405035493133422112284199678361415310914826185898719082055614478428195846280674535073015978829616444783928160484823439627837518823412083904123227193459416159903133796822829758768842795385811216769399996794983491862086448379227447681914328107142619903912754546428365663503288943662286932857090075747789084174101295697803639326361176998613015108226"),
           e_float("0.00013302626963889844717423129807948518006557072628711755725278464192927929152028302463847837899366061550487568568476414348176934001645026354772725859543813918874986991031748398601410786520445811328999157297755062144315435792391686362577021845185481521499609175231729733030529766693781003992275708576061018415023584118995383587400731376272537764699017164495086283786544844147779776869379666690089061153875"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00512_bessel_jv_order1_v(const bool b_write_output)
    {
      return TestCase_case_00512_bessel_jv_order1_v().execute(b_write_output);
    }
  }
}
