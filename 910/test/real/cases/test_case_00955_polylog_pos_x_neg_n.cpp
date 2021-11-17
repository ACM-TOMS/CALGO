
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00955_polylog_pos_x_neg_n : public TestCaseReal
    {
    public:
      TestCase_case_00955_polylog_pos_x_neg_n() { }
      virtual ~TestCase_case_00955_polylog_pos_x_neg_n() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00955_polylog_pos_x_neg_n");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(21u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-(k * 5) - 1, +(20 * k) + ef::euler_gamma());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 21u> a =
        {{
           e_float("3.2292400765214433007016326802958534205819860940281438866625246085348832074107575240170566862854924328940159328216024511589292574814791274917474821597288354644603242581882015068321569009871245280272279619216174050738819437313899886012386338724581823295184640824359064764395582912668144943135757902787249720011402148391811038501739810696042145583938774871153893315402347053486166256245558360147892457475"),
           e_float("-0.31119700596443291673620025994101921969105202073091558727774582884410669457652077287560908603540234197941013466261699444810767748232557290763447393242480408996708337567172015451238221273565495667050868323580005389759444876307366651424812455181026711993670386245943106708429774420214936258901526388690292815162005929237069336036391507006240305961058091466202618593262727964681192553251817144268374389102"),
           e_float("6.0044544135551468950598337928207124842533108289995681761733346031604986673124422672950547721070785703275439121997014477519153141021149360272461021163886366778797612180535791593666575836844914755314153605077252467171526449063673223262043163485014684240330238913406504868481348699571377412521973388092844762036651287596509391625170280844497752367009464691517233581393268525367514321553305879418373149652"),
           e_float("-787.4599279914187890261658911299910684921218254069026707000723564787512889942062486848311302455218962943570652768439702440311832871730879911158077976393146384229365893243926963478112652934313547176195299555180781354511779405666811329497077847947398892632093584928093119733059613466221621334485454574880820611547203352269408375077480840117890359106215944921953092433052669456311354510058053688288172465"),
           e_float("376554.96865781590639569292319915472359828963218112474905364837054890720605476998329697995745669366759304940320079230806320076053141526264678427449581778768061515640603559710129219271518808500531717718341684914083090654504857594530395504812846294893799580412138725853703767615833601361154041673235848075396659950932406271243542820978556885505331180106798390971887048354135500367342151113061384527505088"),
           e_float("-4.8232542237190160978002415823516365829347566933855322165838394992272303664450001317065168750357892994118661045953525753104717130544173094254592366061923131624558278857729285579102114655016092915060532589841557655319014901224814220608343059203163573351540362776145505387324580040606478819325942576234782129573568754348073972989883046258091508208140286509726315330993053113354851298472986924354626799447e8"),
           e_float("1.3729092307249283820818677771139433866271240281342347780965412524133919399054164804250781283080350664645863487222725855908712390877638720015755506508805658925624050183701003946435022663129425229231366656719019872923022849400103107711579603495174801519481186867535388172116856269618950631573044264373889263885814744033370940007756265171601245348294069780801093790003410846169513059872168541569569259802e12"),
           e_float("-7.6546046599788949988668541111219950892118598395655147730655557586880694388169974072352979209778102248436375369315470626258794724033218412941069325124427477384155414602867745901641999123033731874758266541922741940490544659505910802109336856416580499688566702213424551580444423254499278088437977964396378300178850368379886121446164584709353419516947581964864731107372547517025807607770207599885712677803e15"),
           e_float("7.6304884379008440659379836821675178950153142220040806299806127905956792551609818459890049195430910545124492565727504999687344635467389323028016222229602809467037395648805456622516051567836608035650498800128888121101667410720180334445470763188967681131213916663343372471291964804571023622609543192428451517764311046582213686045209489819789345285052143424099797929855789691263675219043839478828130841222e19"),
           e_float("-1.2691047332826260806497535337740700870184538157388641073064858191765163862382515548696269475810259946438224363637098910966509262815317496100956129986930532322719961953218798073454831584736369685233284905339885566928427502243940194504828056435275403035479853362044109130525011594985035246593174203610505526131697202117425319093623407836078231832809743452565678677263311148231021542921275572261218288016e24"),
           e_float("3.3357249836229087484417421623757454843941921448324615115090964723463298164565678478855237968147298650345899589565089089324601167699880879933728546277346170235124487746454335219297176091950692628071247913852479958097211501271368816535401126033798644443249295261708484705020800848710602679190928153755571014022711784987553918507767473115810678139024979021228128145675313807683599611077147447854269171013e28"),
           e_float("-1.3262678222543869458504927632638536091586919350282827131976790350060393995482183079128185640056246215342277311132077429567793963952953222805027435300591478457482829057676048680333996753276409564619515862368023204157146135759616542537681009697683380888890187640722656149805350635465425044256229676906357000112981406271015957074438335051919326737015527446289653070057731698253177146552207227786113967953e33"),
           e_float("7.6943439957538739101895480148455756087883815436172711576992670715433468483761176953209354587350962064411964842756320525451626794573527149596131426641608347697275298737687338659111070973957187991336978171474389455940998265520815146414578373652289964657828328430408996128249384767450643294754559758795109793496115322656406443942062270904407411358976217186029026514676846253527436103753648804466013421464e37"),
           e_float("-6.319667964811411889905356434454616889217933562641145058214643533783500240189836484091471508653187079131571254053042415610509959100626936045124015985253609387135449580330340954444678433516627084266216734795992826432808178882142486729298754869407112706296542531315691375294658407634890715311087551825631421987557708236334791417303226705221447199486199710318639925177928808611194789470222296570362080598e42"),
           e_float("7.1621381866073440778951526395948584098870685303383213280069958366650150801608697087206961394298354558979887936536883435971505897155796999701869378138564553746367776786494208283116462369536636430291109094979057194925905800422554066862417270653573816019872828161220097629514911829519573066896607409568024459655027498965049614987334600288043927190673701781668261074002242648847346458737488952119407949943e47"),
           e_float("-1.095494104677886629882226850692929867889079377219258674025009684048476797878138700807673790253526239404477348363736401709891994696086867618391473359716548035770053961189847360087912945322638631406261012600487839521378639899860979019195032811477886039099357580526457883113118882180358115092850286098062887833049290635606206644923283992025736157791791368733148978630357257989889251675681732298580263261e53"),
           e_float("2.2183796498634011808762950685523123917652262914895467140222274874277691432401693025439921756199504010007175905991020890694228973977906742825738553600162674019265428814069278848078510718330882521286963319224732985621507622054323411461814714242623821740991715257551395654329715657056686789137683136637311122130416634628997864690613805334556022632341829869042033537674347715509988899991474295310217888276e58"),
           e_float("-5.8475967148039904043011575649184290209724968914924686505453657125263330064557978928994375252531880702447500716540169377689434084783363036682564560476103209084181082056705201304871513912024719577293868813123797589520992652289156091062972485514049862730373052437984277473290195815626238269806184550780202927901049375957242684921689496827012911906327512953945267532253874601294728246728083896029826274401e63"),
           e_float("1.9766678873976415545341417895880135812812359207906169401019117840216749116081067903665785633510533459151144913461163515264951953068958374934148808779775397290451459151725188354881371732316679976135890941648910037403031327175225358375343246784102834462651364717628880219212216516630858444539801412072138666182434515142691419209299559893498459734001440226441710688513117389854377294497019948358639457643e69"),
           e_float("-8.4548936239402309805653202222348219346198339677489669376950278934601661771958674591514540367887737679150663210667849944920013611187431928435074871222202908950756637876448805790002512305958953085455729525480537224975981139159700814417456466861664775141893673426003001020866851137649247866248276205310746134913259461013367048655230626391974919640382859163020480329152692873760018069231159516147155621157e74"),
           e_float("4.5216764760906327485790154223904610413143535672019019149006373764231766926506507103103658590172087330073557130581920849827383088916540483360116122306765243587626626177662928395280931145978446199976931271409657713977603590333968560344613799491156295392224794439605216180835920076874806786725335103806896712339820875759549949880304421552298242209408238779351726266983373508294060919134580572195387651841e80"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00955_polylog_pos_x_neg_n(const bool b_write_output)
    {
      return TestCase_case_00955_polylog_pos_x_neg_n().execute(b_write_output);
    }
  }
}
