
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00742_legendre_pv_vary_02 : public TestCaseReal
    {
    public:
      TestCase_case_00742_legendre_pv_vary_02() { }
      virtual ~TestCase_case_00742_legendre_pv_vary_02() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00742_legendre_pv_vary_02");
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
          const e_float v  = -((10 - k) * 13) - sqrt_1_3;
          data[k] = ef::legendre_p(v, x);
        }

      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-0.079031346496558902455028400487692148635948145433573812739295711921530436504368687356685874077141310777140829946981826802229636392722640310977178929601279977055653960521561944017472627589669910285448593850845758946660865271685759381728984476485886486214085785283606780984541991894858918336758037674571175470493504349994958771207706228683823718688928533098513182902625023937121396332915507852287354172081"),
           e_float("0.002892692441277812572980731298847526578797364641323553303943079092448912333633478418055859852884427339369450000702028773443188832067267077610784572916738041643298974380313992533504613248677079232761934769476013004565106864088939887755897552036209223796225934557846299212461704292660298389804587613327092979614492992099109520713935183670391845818019926868544774691629352468234479610658921663806866466268"),
           e_float("-0.079927250605124830678493803423946180389517376449416843205135648642186155396196404605028394995473443143909524372133016704872400957652145287959176086187842777399637082554860342254311089116464288419639016845211975137130383997225653764914931509259691122288521883321109285707349163336568064668126306891187873932889151656249089325979942903174725183130782240912617536510451361378846440620498026374868513901393"),
           e_float("-0.086630502533254707651732450026104710556115216043819949404744394376322315466486350032652564844546247987898876007519349985717300312475434853402397587378600093373422899142409903290897504599247555624299983954787706015316607233101079821793828308400422249163458169824552628828525818352855298091951120465743889541564860409639483172381667638631537987695860362879969661728080305077742192224091825804267228293086"),
           e_float("0.040560296060726995937247013824707124231777344947931740281246302028219858408382923979705755993559044564926428576781787650799788521528044163216351233712551507441422064178720210256365968842322889186648180468547136230245323100274295797836974476173176435624570171111810749217376946995130826470182768961897512363321454637023203883393079959398403865994030680521262377807626937075975136304589846768066301008157"),
           e_float("0.085274177982914442675738531632847912592698789497578322190162549530941184445117428923746632035603879509537141048742660247483649629727903369611844490320698208267054257950555598969556015763095338396862407537894460030892447894546775138015533204654965115399186239401433914607735938684777805463327121158235950285380029751512810483503665167568844046440445577005901442355498879048742017756475753002322295167349"),
           e_float("-0.016861674132557930062862007785485649937220788142818804363807534194392115871673189016996824223882436368239344832763530669141046005571343646521400069391763097435153007914073408581787759871605376098588126648905269648021295888841820258344957805754646617812221319875110494712583007762099858560032615634468383829408518770429220232076674705900439580733965086844756424606259285008198378631962910531306159499137"),
           e_float("0.094809648618166685454853765551270998335997264141583586063254597154275368179861366758088894223390343763334421796973653538382413352608958553055599657229319770986455585630605060945889582992831885333196122586543463973680059711652581101081950799473417363265881284172926773409888817265410331644634289964737222894559165837663714510175513846946166188368066043424968049339430335650752803514403585033280726950139"),
           e_float("0.003165058659104611292760660148279024392448063309913760847507823425552818807834512202709367973918455309958924136276276465227165959537199067575990408718226890081625377745820347653667329172134787957436058086514479447724567624312589357744133110178974749270315809154922582910606674231610589497896088256647386317551638210956876572230525494772515413911669903021279405917211015724274723844130123837543814950628"),
           e_float("0.031455748327928888120446725263392954796302037483390426284190489612629717206615383749308934220161317058232407748447540495121087095010124093673150969888789743558598455125081216247871953849169416981476567779005665293376812624986607242797519454102996710222529926269765588251246267183423490288029077324925069656493992755834635399540835958666494544141971342170816814754746239081678363352009565277109595377579"),
           e_float("1.0007053945334960366255749083305905652993693614001548027423116259096385859296367372874412211983564861520849560998619397586530545852475248171737023144061034406504056283811182460543648403833072498205264822064177195489618645497996999370920755549793806991088637232221244798018167402046648854302501459201324519316809968628951618976771264310817239814145331322242393294293588148133700889016608020153606105609"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00742_legendre_pv_vary_02(const bool b_write_output)
    {
      return TestCase_case_00742_legendre_pv_vary_02().execute(b_write_output);
    }
  }
}
