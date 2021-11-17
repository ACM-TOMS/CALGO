
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00703_legendre_pnm_vary_03 : public TestCaseReal
    {
    public:
      TestCase_case_00703_legendre_pnm_vary_03() { }
      virtual ~TestCase_case_00703_legendre_pnm_vary_03() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00703_legendre_pnm_vary_03");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const INT32   n  = +((10 - k) * 13);
          const INT32   m  = -((10 - k) * 11);
          data[k] = ef::legendre_p(n, m, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("9.8260088427100394312036519470721803590335581176190487053061891834379319714857160044685045034656728782259534672011845365666789136022179636108661786010288733374899668780598401735770337410068176984432564866003460736364967456420569916650395680750359192843902021100357204526385464522661052879282660071169690505063867537678932375048009532774463084625737680205533967230111377363516166026831442138726316613657e-319"),
           e_float("6.3320054137037191690694045981950669334739239916465393375955906520762140438747279744933433066039710003622213295425110708541724828966621451442909304836501119426573770020855231193234825403192482239720789503760599730982207174828792846119351338593824429786568250528424558199002698020523127112033276223121513694152791116462645759085177686173265084593236762099611690601852643941108353842171334196034012929844e-210"),
           e_float("4.4720866770525904724194504817222883455970188862024892743595100132192922950817936365220137666362113493454021240501886244336324685726589380294731735012023886434368989382489290581360483933260112466886609170808127992839490172969952966790204647049995344908189576439903143566599877242850833741209875073893720883624113547233241471031204514490724321421896115968796799688382666551964625517974537117521573269563e-174"),
           e_float("-2.398489108606685134013146645247507580569340718969692323636614970203748511797313982382647705188157900232997957389675078681526887228412736145963588838020172756089792150600007992048318340425148331653651287082924373735280051826463111890410796005938437400686345101890761326254354262669771458700803637156272522400487235530017759978399583293940221335824685736876819482590925936672916216455541664872151127318e-147"),
           e_float("-1.1912609411581729651084352251696228218435406553176442354311907791720184274444915845862042137097091130806957650248227755938377680525828868261169670408922755485804948467665071570689055804385322797388524287012211113209534606060442112677263957756504991564525744456988327868794584920877550793270382429478289762411844896905534659982352745166201352836294805715064848331539926515064741177009106698812226992657e-122"),
           e_float("-9.6436193595822303787829815704076173743024566433416327934344164108825202966730610871634070465002750570080518226842512507181266838601643231027569637609984155707668622040097154307480474301346071107644736183451960055512528350844436082689225701336728986359276858934771825164364103200008403553047604205068203714929862990140667998945398482618393208027467442495198088302263439622209378716562191562601503592017e-98"),
           e_float("2.2398094997975152300401521251169146301207737150696746342582495379174753408842197879041676988435957113800149880086950941346736960812072470405686379074211914844354340242717720606363229510618195384031111393383023579831457343068131013016674396284245613386436948083938526386021224787778501954750323158898880265771695910153649669506377196387475670980664584435371869416497906143148061376698277373987041719565e-74"),
           e_float("2.2760530572690410791130910581290953876545512322813171115671677876389306995087247829752363873217959193588358084029005642224437152438864986742275821760580531602786403764823584466256478945177740639802759193098488199258186325744323462773875840267836771457690962701501270755750123505753849213566938832212950173863380322788683542148027807247559263578889556908586960303081811249164250039195525897811351621967e-52"),
           e_float("1.6809239017581002820567829791584937494213013423296013881489376582431404827124252898589141946471363026396138256257803962089945944465738809507701050899190815117700038335796525686030077999111967848530411174212620932625988540312970935910371010984624980765651577438520421625886463851390552955976428426126958282079329763170155357037101605741963144787075651575123815783156523297328392813000891712768235404208e-31"),
           e_float("3.1378712655453888131053437903146423878251630261701479498401470085524256333490442496123531617974913143103938095653703789995095743973188858733806871736515896493364919798267908483813496713026559166957734727793963328132663426387029923291627660339971438985466716493349357185643401575102598066803086956310858987220955351653509810496390161242439941118704165940783872339378038102467166304512343667263158474848e-14"),
           e_float("1."),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00703_legendre_pnm_vary_03(const bool b_write_output)
    {
      return TestCase_case_00703_legendre_pnm_vary_03().execute(b_write_output);
    }
  }
}
