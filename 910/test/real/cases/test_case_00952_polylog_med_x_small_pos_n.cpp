
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00952_polylog_med_x_small_pos_n : public TestCaseReal
    {
    public:
      TestCase_case_00952_polylog_med_x_small_pos_n() { }
      virtual ~TestCase_case_00952_polylog_med_x_small_pos_n() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00952_polylog_med_x_small_pos_n");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(k, (ef::euler_gamma() + k) / 13);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("0.04646427478183840308813921649717201402013441152388293171046136867119904916935799006845937723601732929494021991180840342116144599052390837522463639085323401012208302957299656618398671604731067971009443205896195766206336034406619815568163949199671921221229076670180093610347232258448062339113524243723582869743763024506236434011029862812073982713463606039256325256796219967163326537758904288405664985809"),
           e_float("0.12933937077078092967083526033439449988811501270822050194784014640791180110366694123954958448404520372817151766464765533596600083426253658717274634090340201757863696480803303192811090314999413994187080925274948362698533307970227764786307524319134032406636744484566895646018370290220115868695160319384018982695908275882357963244793618979660639241869394161926045266820780166442419596019979295340677687366"),
           e_float("0.20904935279470170500031159479636160846734507496400570274128880801020101972918852675671880281475419877243375133349678157509935136523994708064188845826043364189826695772255626260036921405719115070161391011371966094915750012399722988724133648926687772348179645260854737294017003365631585219604621281968114900001870829569454774358944314476051941543461209559441256643409753956629389215028229726405269838709"),
           e_float("0.28551161180395556302417026423908631903283131314069304391757503343440673063003697701579463599549600459270238250077091453777614188164596384164070656761489455014944476309909829052178850890921317547805594543592938929459716031475726695289904569310302959705765119025310421428504188436361846590087468093963915300579343465932122036050399677821491971623860161618240926751326578420763449864323057448251813097507"),
           e_float("0.3604510188739002818629455321242331421933727875450463654924765923460909948070723235505369644960809314415459802413180771100330243040851064279989439088751758793039636496324806534883530754899299235498661019543099886152014690224528139403156509031622841078232768679385650907406334791284688210808344991465615335449754460254999041268902904611581675856645731519533853985480687044974340177869947238458554772428"),
           e_float("0.43513200558834796553841046840105483620153198593160894539850243792138832989448063016037259738937623999508433864192888374760173288948552648900171952010066004205931224648897881272401701943283693169408993431961754091971528836700092266745681199933940266081605294730304416236698346403324773803019552055217863481722689014014760486505145808775076362630529961646028237758365764226894850885147796307840896509752"),
           e_float("0.5101354986261104127176530112914478115461436444842177741226131109455086401468143700646016732120477859126118545179285719062495746328919866188651191259779108070655666706655778383592195797213375458807937781539881350954489236190596060000787837956596358689022860536260877635405883749241515080435009620660057555852941006327384402419028077415722962259018819164948370691578407410120386630321451180693151720011"),
           e_float("0.58561550003989085726564693275426832929681830792236087212666319314088450719229598591368311027075236601967816849941818961157968115175912754056702686023295431380272372672905264473910986705788051641995634332513667344930082413509408294222937956989086700132834019244918759029399022505793128313738063616869400342776317828370235930611578649480971587850812544108225438730090838002149770289799157416731461751374"),
           e_float("0.66153332808815851897667610217650015731061089951058011462973361590818231655513152572938920275971059363515346748489968926335143400782341223924357389253009816932602219571877503232647377906617376175654122015937384372690247014938015926769285517927309691239573833483055693553994980876723527908296864097513312981074200990496497689646455405210080308771652739100664362183937316763974127505171681553891290849191"),
           e_float("0.7377905048641279800033300192195675231471349690052736587461519646215755598290600473944317581044075819771738264503062947213878804235960792148741526442986729497899383570902274851973964235813070568471658055168807301883261176278628686604158259310371933243605426812290938071627884780657871514853354379923019907416532509262073724549445673213019654595276967406128068048330970607187953859595770831075644650353"),
           e_float("0.81428803759446583518153261879450373332222005579373905864405796690420563617592471007733974201219113875952868942402205626043844986216194784318313694415032197210142999707886522948229210454379366066415640313998519274935340781144543394657900302967078890444846323473186801849004055574817497145756804891948165490646152935144750583317657688348736893582804641687788643976367074510441063909777775019731143451454"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00952_polylog_med_x_small_pos_n(const bool b_write_output)
    {
      return TestCase_case_00952_polylog_med_x_small_pos_n().execute(b_write_output);
    }
  }
}
