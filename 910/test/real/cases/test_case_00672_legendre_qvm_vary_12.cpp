
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00672_legendre_qvm_vary_12 : public TestCaseReal
    {
    public:
      TestCase_case_00672_legendre_qvm_vary_12() { }
      virtual ~TestCase_case_00672_legendre_qvm_vary_12() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00672_legendre_qvm_vary_12");
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
          const e_float v  = -127 - sqrt_1_3;
          const INT32   m  = +17;
          data[k] = ef::legendre_q(v, m, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-2.9500656758589144493762001663582785503353546441927242714259498567623444066477098347647365005914293094811841611098338351577112635502815989073960706338271871228050481267137695201210220562338267132103374779517815427208811068348416325563727678592411693142199544517111400522871238589194416577661692997447807766832526898327171597946295620309379979274570591855598708801846275421244499419462568914266650911575e35"),
           e_float("-7.0877037578311796684074372403056820094667935085939345492248254423811415029593251471044437639656120077316882761164761028861941809139448373842526478150007210451142856150453693179902972947592703411480375176171018513022254196389953348480686802808299975903639087092558447691500294085621521331474124042376483334290257652808006454820272213783982711822429647077035683637220497849129557391222082694077650622675e34"),
           e_float("2.6844606514618102217806874349816352152058002171774702885505288254183366606302946538345972912894470605313003729783765181513390209102113037265134642406289420689338485479986131477141474203116322912111683520086567552201316080191973465500153779662091187155080075949134051925607089594635150750759420066089922597923430804191154852972122446661412189645681660534163682162842718903965695848687373906759642841035e34"),
           e_float("2.7237650837386448607510858941795588619265425248911326326162338652301837547157050331161842447459416077589541651304383022354044229625615701955847091719846049551031539688312590906481864650136801566761058963677688816023122155139038457659192473898889787394690527293490927922400497206977990257485059114649563357436552447273219598890211014375581929855208233942866322664558686351330725390703408581125992564774e34"),
           e_float("7.0124207002195347712274937889660437171695802231918825098501955990432062598324290423534020911227511578177562966798469886070243674156903672342441653848183926965062208104873866958935699715081964120470659384083297435849128699717136015721322378573718474433117421031322864148764169982395155890261539370393701779438049220631204004004803444396871508486270964474983879869071996664089072092023966262732470954687e34"),
           e_float("6.6905698143854792409127199412498991856341494686999124732651942768585949110049706606962485881840205474618340890630520926469856831029115168305477163116535659767673085695713020339045509828143168832426055009243037755343019907386577910285829568454911571654513024688030071375684575573328470203107311063414252857013126317460710518982390489960107908139593652841113906081384984661510501763248006129802665992218e34"),
           e_float("4.0006287551606738279137855591684752717946170647909454795943267145863852661329998522043709340793563866772225630974797886236019760152511479889765798180853064047256947229192355872741913108282286119073296770037083678962652388268657258713295266816131847510881273576015784587427510112927873337585289550283358031908834029666872406974507142499114667122761090129419356211729768696064271310991294328399094549198e34"),
           e_float("6.9578951331791310423262218718326385292104236158091822068401862597644339489693617742165918534569073857418300694565301438803424133056799871475252356124981379990989288276067564972819369862204152352692560546873112687069999010785104573832678432525317455691581095730926854614335915495308656750575955604995318705575588628622144499962531967083534292764317945017681099119119887916021087545690992417705042478557e34"),
           e_float("-4.0084341705518916200956706919119205136812712940852009468687721490806369970729554842829210696374511809196880103148545451520164949393353078540145030465649517697601823141192749192353685231350224259579980359151283739418346081244736335521692557351152519072791742507907738044999623266410340381618513706674875870319331717441322279266954796055814130837607018726474451710071140537363579009430520757655363706332e34"),
           e_float("-1.1876331068304156902018838366414984234479218211646088993530733856107037343456146881471588816277174868742381439015153877709814637361606578028951315905066085773648789039586292068654909275950997118262178255997150437486349840026706613872723991963162307309795199072215668401852176567278796222927297365094595702363248779777473408274596970435590418233831394966456423707824865613590925173447593769251378065673e33"),
           e_float("-1.137010862713865624101815655201007149924181137005721128852789592450800956662147231582497650262457819273178696977520600526800802236198946468336008387205827907783610369438228997733673393459928893455178469077347317452941873194433781303032885368163459247825338300796741999293655834202535441600571425031228487901203077947630084814272429939300832375204703883217528426963867367564690819887685591328231228257e36"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00672_legendre_qvm_vary_12(const bool b_write_output)
    {
      return TestCase_case_00672_legendre_qvm_vary_12().execute(b_write_output);
    }
  }
}
