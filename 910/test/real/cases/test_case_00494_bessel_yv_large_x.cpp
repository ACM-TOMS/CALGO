
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00494_bessel_yv_large_x : public TestCaseReal
    {
    public:
      TestCase_case_00494_bessel_yv_large_x() { }
      virtual ~TestCase_case_00494_bessel_yv_large_x() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00494_bessel_yv_large_x");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_y(100 + ef::third(), ef::euler_gamma() + ef::million() * (k + 1));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-0.00070531543071161967511703990622030559951279251871130305529403867138432923649969827516939748507438798685392599893752431008334302273104169385942943742391334042701561080945503243286426829802742588153828881221517035927271326629620955289813961734646033032827289229104142398999987904838673298394568560709487691912023940630915304684673000555078516133009210980535036021328499209051525898337811265591662643916342"),
           e_float("-0.00055968824518825180528493483777952647467781554625980088250251812273867058366785842098482776757626901054577178816104429775844235324664041604950166415940355349627314856890335090813281269982781363553810303341345252432556656176414372701690712933782460444799317925656638650635595390461797208534082087661948752243745441904984489637722813341292465448327631128598059338621937473523059649969511564131754129618478"),
           e_float("-0.00044831711007771766792537731838630287103880736934181308931789077831037888931536978262569617588923344494335787397060501130831375752190922506244185609638292840498776210937807770973443103630749556852310939743985711074996199130636946611928316348213154932549207324707494710923685100446569036917675174657302905420857342768757957821619775174823800365614846687977735680169781693913457054112913090400867954611571"),
           e_float("-0.00033150105218306752917814412454826460941930243223050772916106163833258903578893208969970348653053323853219443794778417108548864121660647024670794590133645011104482229415588628227728524538242217933135180277938349948195333992736716925523536855497544075769100020932596437631724117421370437005491347391153554709785581821875752285673006889660785480863550163006344173686325941626591135552008468731942232292532"),
           e_float("-0.00020819724592377249396124961594685780994200838492860622109255308904180125210497296520750291195471596482513863884251492722818452915106083721621131319609021233675767359859151487927022106196487739229245412670313945047541741497079942565010058672906348640240149834120488941015445645017194531814287243292047780374739305395963198673018014436004686452787131085257466955569902691693158430064119128710876568244862"),
           e_float("-0.000085396396920822832076606703753310261894319102491838828599241311670134936067830412334372834474965307588638542623131590768215099076857376568017757659617550853982745982184429474568756543079420547872530911046017894648002007678200413272247859944963814836712138921807148757864367486187116523206708366976065755562237808807366375898164293412843437602850523155987581360666263925643415782244902600639037779504181"),
           e_float("0.000027831250145598544034081566721342421066507166391389372446491927015221806214855910477758543671570320145181539968463058933349074197836884789444474860685080876130855831427274713571311859208458694301623110235168928967549454595096254069715495610504075939611431662647618232381752524152834600912479112004662614076619059512448814548106746733869061543642058204664721207544005914510451962652401086847796934024226"),
           e_float("0.00012271999886840282513971411402513818094371039291640476694484082843480937447962265445670381665580696736726274115109527033960808061357221642607430471554607977861373343438499187412732479514510757741365268151720866917847350741294153852199722876795174430439530641429505858269335149760977261576812730455897441939075409194326132376422801407614859028623916875209812514783632992795926392631709769260503555530206"),
           e_float("0.00019221149195460267069784827002338867901485711129610732824399665043839463604568059810166606606321675085893605461482826845056036651010360932507546585136075050773081838747510755545120735257126145084225990963317118137737947358058835714998204044886881828488284203564350909900549085165892503054444856165287657364571899564940633349711649531728577652040411393886768264199506330357236971667608496436290569128239"),
           e_float("0.00023185502056287209720826407543648364017895317162411443762980422603852707859156841176660293942368460704263818453610162003314296328551597737456260217875692481578963517507167205144597092289729188758196566257348931422944159092372837364297836457524031627070260162376913042860511066607551315277339023038103789708815081951654000490149601013749121689527234542441660223640786623615238406673262815375853054643794"),
           e_float("0.00024029577571008146108854217072749707278926736892239984822841523306847157725201850019917039066882058561276610275758224466922249803951863141655202395620821982684358400407717838421155086006429703565096602451721587092923377083444792591796484457564869878158413804710645230531990172029458510454594035296093178637905366689732890949351867402824494277847072868682220700038417712341197146687036796230626844312382"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00494_bessel_yv_large_x(const bool b_write_output)
    {
      return TestCase_case_00494_bessel_yv_large_x().execute(b_write_output);
    }
  }
}
