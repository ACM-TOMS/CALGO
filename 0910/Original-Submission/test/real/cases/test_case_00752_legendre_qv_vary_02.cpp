
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00752_legendre_qv_vary_02 : public TestCaseReal
    {
    public:
      TestCase_case_00752_legendre_qv_vary_02() { }
      virtual ~TestCase_case_00752_legendre_qv_vary_02() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00752_legendre_qv_vary_02");
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
          data[k] = ef::legendre_q(v, x);
        }

      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("0.25001253162842850821457968709902891548934693126925373790243437880647698288712662926274333957673341275843316259711150817930071056300535859357524465004593208847885912790780573727569863138074061947456940042309909545403840491756269044898806515087760828557986596574246587382392688379326387038337701867489105572516578711116508312663546390634665455605004111672814837868585855871084669948359351240765135838938"),
           e_float("-0.14627193969836270234868217497071831743287959625031857696170264125848001174730960022351441188267506729919301983481226284481005233392328571460801096197837572904325820907819997281925097297103866822667079539320750836338643444870424760954986329249470804720095954654149024768269417671766093816708067666525980483731459930018456989881466511442082691449484513678345912020474290774961900160447006899791815536098"),
           e_float("-0.1170379665878396442473209316412040930173105662028004840161904356418243458540668267319541189238480268700043573491237652749405414278855717844807807661219707540153282929375782115494269004603789835541176794307507412761001378021984385068354957102557619539137511488266831847671926547351231972837874072293902836007204476933609217105141245386906996113647012250985896748865975463378961756038478247503615889079"),
           e_float("-0.05169164385847892552709965602352446214218018602387494559659477054142038434214813797884568972144412193462410922856622260721320431667512518766619344769660760546185139530920335256637973854080631184541389301056888197170954087540242762849581055564401377435869861725090293011622820674081076367645721460193083109173797874241333162257039491063584494507274641162314183927080805263101459337681760091231925832086"),
           e_float("0.15984376277214734580451233692263258448958114396150862713552354030461007281985525639625757983801793238066754526936049813662285050946975111663087562451590108290373412229830444324972121203714540602280364442804021772905342184385738796640515246834616623296994408693573263778579026391959205363992581230472672026134238085416529657979557274441213816966068167895302130386606498596104387706588955972864302265679"),
           e_float("-0.012295563007526056350935075114505620906572247758847957272546090827022616896134517224299527679334623743459253412509770264480857200107157153338520991217311275359789807770733261231101151379721292694452600579314709458865128388957463038502310744245130775976825875818844336163478832350078744850827357004268055993671278432113702737304650343542152701892160493510356296040984313123992917902627033735044904161819"),
           e_float("-0.18646783969172457599612107087068506790424998240452468789564655325975638600920101233571455463131412111913372910392888452357527044768622902740018653629879853635250596431338577449472972513871622425800015861024331033338177907041753828239310140052493423508700788294467730046553327238228458474954297264044626133178247504078067561236207409108773409906114634863552392382745569140032893434283132241472533063745"),
           e_float("-0.072981411698573197112224303839985697559043203803081097502360826040231531588958848827745352641239244274245399057060100441809103246838439909725117077531241886334199320518931975559553287779003245271660408150820090601114356214282963313805646860675221357756435586407490681435661756044724671211941044224751138251940740398340331829953953691547883006167363885046520298457757502296011541930507812867963304369816"),
           e_float("0.27604650603744896001719124809832289950090141018722746389397693564545479911658187414197666307259516589206823470855235825722748316818335616122525277707445966139947953981437254693219800177033919744802177896157952238763628571252711755394885758888459247385813548357332292799249618325616913963014295934634313765090622465756054148187470529243344519022002811463491373082505915112219590775093831905044554940952"),
           e_float("-0.41692849869118191202877763895423338346849352747559893131750482594628867319174743924981877426282069847409650931523934131835506938995611797050136414928331337985430049779560795543826640404042054050197747337560466327398290638918812420743768169759066465221971126666406468386537040242555478192178053405382962961530641146524079322470062290096072906948633776659260788184219209526911322153296074948793510106215"),
           e_float("4.7537600722871229706417355380759147687043969895436784293354171540929477518750349826283984615712697728477308621140381216271048827692623805450495740005279872788866402377977971212016388650670461543380532966027566959875675892694091104633192546606845191224925652698000225711499023835090084199943039937411958882093668259880154373033256794455125285162363514050672636940561305999823202527283549502484461485745"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00752_legendre_qv_vary_02(const bool b_write_output)
    {
      return TestCase_case_00752_legendre_qv_vary_02().execute(b_write_output);
    }
  }
}