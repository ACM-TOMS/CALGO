
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00112_cos_x_near_pi_half : public TestCaseReal
    {
    public:
      TestCase_case_00112_cos_x_near_pi_half() { }
      virtual ~TestCase_case_00112_cos_x_near_pi_half() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00112_cos_x_near_pi_half");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(51u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::cos(ef::pi_half() - ((ef::euler_gamma() + k) / 523));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 51u> a =
        {{
           e_float("0.001103662615143147017601393734232421549156847219973314118949247990939859507971857976111288683988415478257878902917354105871236439563170045857810738102729287303674642020538722348587741395785016549139854168605590336633274890874911975139498047488076101124170425114828336211817830671418321600962996502656221412102902998671173649312277667706874059520046508702908872747388744526457997404065469953155187340837"),
           e_float("0.0030157041572432190019717764776322147049224887913331212054173669461673516770769148735558540290915754514269097799280330158745508542936966046640795712691698490758472574837678659161855207596563250546269614348515005196657786359646020537231301121492728577410212550350070502426482802490452572415879214994201864402901595309726458026672888041860289772080407882673424767232857435141510615538430445535588448032261"),
           e_float("0.0049277346741750694119109616231402901360500185022482787362616412856233338679208830974124749830793778323997684640684587346403570647040211088376986993615001420565608501262335822815041623056442245053684626079448221530207757672939048030945164948521320401574074802966260975904478775314730002799895049652983052125151436511589335123254582883803928373794668948878049089088843407672795237410750365571990834236348"),
           e_float("0.0068397471757112108044544339281399968563253858411528016363900245610068161564091499398024110635089866212789583403491711476788985949178019078796740862827299289838108117219486357027460681238178010067475221327133396780718266880540272936860185783976937569903778160170225166328920201699827472893824363838439627437187765159770262282930293057633886299749088561182830982345936603553661004526800797082995898813903"),
           e_float("0.008751734671690018552720976028788304981263251327799787402847843531098338033456389098956160956160597460869666211569183483587979772678557230743050583975270977754518475970489741156524119378847305446845331098526962579804361960647344836744631465315238254178161716870409320747176674142306881908614436607977326507761331050644035013756409444911127742401762647885817361924274132830920030936589562922803846509607"),
           e_float("0.010663690172041286306192948782895990881227033104763005629047563262531418373465237876680529005058533981010793759446966868628480908709812402179047774888094761104238131843239037167212206469894485599085717216855637620318723457013401360290343196462723174177342294339229792213396006286331141933952648238459364417067676760977202812721209681621184724443514062993874571816382775241753902032787099435678672474841"),
           e_float("0.012575606686811781116778815864446643746718166744927923785139233710240017201737625088406861185781739628293213436323314376708933534933748436720721568462267317026688561438144113855290140925120435317506169856461386623828941531577425811393721991195661429557079268566156800689149725770186985872449679872155877537432897653904213375816262709009448220490217917537801655330411526526290708967440087145630462068976"),
           e_float("0.014487477226190798137230457835664944872810802211308342685055637400305474021552028493778559102552025380880477103065187688337931713540557562373991431880251553620947836940726544553647739994157300967262941976933375115946629204615918888689039777423610894099269527633757403818072000860056651445780074486612843419789636984117740884711882675016016919064788076097669389652323393864660933705259163296837955978748"),
           e_float("0.016399294800535714798489384358524900105286836154421890562438917862501411947245037986455919696365301666112414890097246479762149456962809638588245466988584093380263545543023411087203351948681559905336010066900667101787138862856447778428938265451537422367650992162004154495705878092019403922006152498705687706013973674800621422751806347308220226407830977673540994341535381202679482907460136241543505771361"),
           e_float("0.018311052420397544372537858201001145083914223029622271478159852925692124442640247100336262041440360412392318902558295108452513123454103504516607498763863728284797202978975409264141138861738184468257576157204552383675442780049333490201800526935793000727771598103735596829143729912656158790924063286176336412427525710457367751437348230246446466833498557166002214984255422115833828533909144072330080184872"),
           e_float("0.020222743096546488827333191237246591760639730704704148795009602627175281013458865557238240830203973690492882935553111519759474835556686812553430134442713057622274348510526931333124437928482001301214776111358848026164688447150395879262705323533149853785508106018564787140458587165768350143018847840661371060583295055600185059031740386698369796694993301139380463919242175430688496936048982929502368583755"),
           e_float("0.022134359839997490880406060727014362125793630864951833270601574861700025448109186301755833862830343973100181180995751924514611578149687732808840212399134735559355116816700346196685057916461444889179117142807858883574426062415931319905031345377305247124198345242221848293630112578508459623098934259229893301963120289071135120089672338617200639338087265986570094025605269870532431816151660553037324591448"),
           e_float("0.024045895662035785157706623778569329505432699504171659064296747761761946910295074291838826981712792368282349553905175996858645967280135417577488810651594776626562252686758463813820791031580579721782251705600570838044158357396828391251385978824608388341556751951593319015721256942139947235672567093309227673005475237474351604086627391812539914665025071874865132198040631867741985556198082915821590691115"),
           e_float("0.025957343574242448364285479040344211549214421772336408904791690210449825580868841229728937909141136116652894475476991471494524501035809956288247010861464302066154254736425460468024350571584780489213056428098560686376999562953811729945811217829999284678155635497607933780499844559793793320393174691097636750027421158473345884329284455326733063073401001936304522649050748641723222539485352644237494855455"),
           e_float("0.027868696588519948373400137317728007608587427701668964330844960393232962659347088604492706251530487148805683505883152882204177723404307465324383408418716579280515871949790204586171885981761699066301009410176801635003961167315448162116236938093263170474354430610732127054993766008385425707409224417795873645410443727472246508894027549288331903202127478347028519320691939639175445464091747341289465807076"),
           e_float("0.029779947717117692140641616959423088975059158964671307393654571763491553590582559015283533532986018614029995854278748039155177088839414450454409571267110720713374288363670804692289394770804051182510495900322915294613479146430103948686457060129295537542317247965142780968717375277053102667416124807987527429828758851341378395347636356414345553353015296944031211364374602764176767114282660874356057979248"),
           e_float("0.03169108997265757234968007548831539092664261419258315028339932724609641378752993073240305448323210275610250474048387518003293691429670190316839718892647475151605557113852966925288722557261748981001717642661075734993748510506277165730472207647763678204813714194082852502621623007244214605646094787730992510318924772328958129857640953691466395972638014259980360897599169526022683004620965496473543834477"),
           e_float("0.033602116368159512696233026049065446239504816931754722314316906991719894135599777826750609185018647633748023602662568553364455253741994701511752478506175671449029859193583332380924298164350627023797218098539441228986492601142199774579432202424360114528105854714450548265654769162069979862894518829365029262734944196495563243510942310738160350648115273816371164862001848102516957052031698646996790510953"),
           e_float("0.035513019917067011716864665791693591701266221069066142988726516832036728613300457006403299601385922973285573124147669717618486570986860008564309641657943172757912451506572783449830586163814477502720341010870119231030511664725149033496023356743601581035789746158436049009960920433265660001864258067191594425969638702851317416487965304378219676062220223886417078628342731900948906319081697480745936074075"),
           e_float("0.037423793633272685069230163289099410548201460485083250584815890545695501760634737348027106921155349564805321222282824709072582234185849564440670000847571561359324629877368620222775856068569278590582843550096096627510109210574429458414697641703261368397349100700647508781941466979105505468662952228018391181441110331804231754679360190526351574820397399161022517352576803713854137820941977691824804918229"),
           e_float("0.039334430531143806170384413477273936914540782490679532559795199052084263353981243128408259086619593938946158294175465792633584421592350011455198124846778194417281005620258982615227234871780679327289232072858751285307955384102986030326321993554199425677956531006546270705693578070953933490084008156062047268010604414954653769711129551190023978078265436680151306319349753427763894284339669614362603225399"),
           e_float("0.041244923625547845099780771389235760645292527949481044094490883777921290140650504322516472637184803214130714431279474404434513167088543202171922178874118967991301744186858054508473413647148744260057423931033177973119110921967581379650616945061418716237102199971097086754171168500651650526166023618723542873290120298640370217463634503988257867984189842393302366909916756244585593509804376197012137137735"),
           e_float("0.043155265931878005673591620105503250213648422568500347993529171942176855000247082147052189275394127432389966427668413127665366489709229602926384588361071126297040547180622989395253518259197742053810419372452372506534313264090033740095357809571640645325965973975968081307185286025492280541542290934505383760407808393490207181258819590768103829347539943429968533891105456889889207923129412409325200337396"),
           e_float("0.045065450466078760596989313842010956972889269629840996618541471076456145341576346037140737554087089598232657340541865108389258265511462307708855932308284610805209084413460614598570062354088931179587984702006847566162492414964136182928014974554852612762940578089297583097434972291633301329615187665009664697057933900728259220674217919826618275818024987757805509580781060051826778978395144370327934582031"),
           e_float("0.046975470244671384601033063916347240271496301008030038885300525267985292259051299166038486168462131096926077920611254966322965395673221324506213593415263189667393727124978454125104192236004303056348320706091018781990803164595039969574510002737990997634364247271320535459089400357647302806652831431190979518269432197656547265417669309719787317687585776029758282182901563328846367691507505917911626458378"),
           e_float("0.048885318284779485470814703449342779874157141188272539012461507790598796169844054420701400747288103863351008593834375336008952419161297787680355981672505884405645257094616683497040471721980241679616912483788844007649447969703512155758482444019770785766752697959554691538291379700303132353876869628124563105576278585509058450081192711918856522653293896204836307264696289474219755720645459588331890558889"),
           e_float("0.050794987604154532871523976044438688478988546522275702591135165295075684834775069782087540212155311599924632216878717363279038406203747968637573994688146329296755023729170505930935028107359530042909822369245609091696704070694133386183688825674952803613566914142157159664703170379754943610785780639205752343103591085978554489192416588288410803225066825839880971435465232206429722600557780654218330146221"),
           e_float("0.052704471221201384879102044313723519214541830844294339002006150295000261803629777930668202713875629368914012613840180174156598399096984605842000917709330893871365961242812860633022132663435078869433262371062200555333508565564704972551921154864927879748108180002977668648111714062992755985320189721615088748697314805428967770533095416911693326757965362169171309026356339280576158225848771071227224373815"),
           e_float("0.054613762155003812122160305957618147686113745063989650607179163157398820817715079449268802189362998719256947035541552524839497119721884536865240705017403409053301625152685151644761334409897199810093505097528650897729279062117617766845257388412635800174656504897022571125534403010151970450307624961920347280823962725005822339087173498538012712497715370388872960900453285934131578336748483622553313648565"),
           e_float("0.056522853425350019441850338124226420404372628229474954730316461224029279793449121148941830155705856364290655296696008865349058650396297667951745639266344993911111272076019522476315254920931469406381088740392436352997178777214659473302034620137440843681809428909274912969436139312159105442160243076856936138443854294531618433726198293972487505166088883354213979143973990964298454127970039847456893860203"),
           e_float("0.058431738052758164976379864942430396555977098734986011309321413104802080798824039565461521608355122102022725461220269967251619843394860641070621951433765763369761658916816880868188275032453137184211673934079277285593935159478757148621410070196500221526968060556586710920750728824112282461003224923995158937905479480579208549847478185757177039153501786271361921802927586114942206193581554801605680183774"),
           e_float("0.060340409058501876576879058406717550375449897922866349024612531373677724965692896763451690034750214807393086060232183962741893455630800282528706094634312497867863194252914455062127544671909424505367559503649890477570116767021509485710626043408265846213906257952113003324887350087083283630715857354116849013376841533857321740197047217435805280085872071354373524395263864933170034771001607559366437835689"),
           e_float("0.06224885946463576546133123915706193746686049729742751789281333411875145041642253752249398319274836506892916659694331849766284472811963607782488993959669913942440995816339529777656740297867471191540909520680122909159182411803441391065468297111274482089347737322313334144161837035757661055283810143238846811496286086926075008068489724098160672664550444107807861952850135669761489736449089327971794824585"),
           e_float("0.064157082294020937013292141108273585068933185681722806586173136512969086561151000586083228827704474922133659368208499428517560937267711287558275856325851735527047525549135602809953924346232393575354505341516408877034795144989051337120442048812961016500803766780090428411925756071299228254760633032932651170483895163866111455935100454818941073761804284777739005007991141161823471943803047351369535073521"),
           e_float("0.066065070570350498632132342303164932774411476452581548988386145454309844170784631219900687545375669726062411278210148833761887075889382631711319762066717518842787273543143969874177059386003262516941215168234997245420735913778365997498382320917159517052861653706668051150014889560865161606396037665048820655743861710636504450723955435275019348223982732095398949832565423361714391500479155224225613768011"),
           e_float("0.067972817318175064541548243754498015429928588516517183822313264446879635682093113384264993280408482780193047687370305892309106748875519486540879020142902362307237689684697063792225838111742241714789823601251998299331254097504567214769822096812914869886418586879378868131071702869283011709161582762527864014863950459121468893761825841337802563429112718615110803469902519051203360473993644451463303576859"),
           e_float("0.069880315562928257463098098362562671709183706398122146969098129179472379455539984270317889382725678551314678030902515782188203802362460214070192378666311369730755912670660970378408284876465192356340335277225457776091515934608751710356445501169982098020701963512215973413220903547301855847279762573955110834619967590588317288171580291576756417265600792610005176888567344924086299458709493921684370833814"),
           e_float("0.071787558330952207061531053207739391124226218288071373402959685303573894071937999096931283415397385296306228484128726125223568292760110763075739612174182840148193517540379302896938375357744754831885570482299802985278525692833941172622048860945698352794727105726197891701070706523197225463566498764869186174536278875440544967870368440256565084805804998877012238900798874720340885805165811704445925369809"),
           e_float("0.07369453864952304506868897057861922170414313661538863465614770246556638808351254050794164760581857499175318127117745597595046399453737415234966090460819383278338184265518294393589316067434296307320873514873375534499000895507092920123162793722742278902956263178675873980121089525814840115300730124457959632013360479776554262210801798434607194598511786777350540581826836696200671400628870250659962729555"),
           e_float("0.075601249546876396992772935963040077834759122887264324876826746597298507042570795315940378957369646666545616927930046677179646722187557860442506025110541968257406052679604681892086517951648112201767029044637142334937786145563391280677987040697732028896674786378870229133792190093646530752257521390201161387698590619925316803401938362974478273195688740273778224688121208296858072701246429755713258468481"),
           e_float("0.077507684052232870319778844857360483087988665809157375994805235390695227945273826842335916971262245322463861698987988491433092759415620685754813456247989287157909601376742141647886684380473379403067343977633347655683745871434202230624968375294472571053679501320680048739346844244553807193670436672857163839591129928936603341033670518982144557944158002407590669543043048232580344768659176122645763500513"),
           e_float("0.079413835195823539113919284592789131508074280056913456617148863190077260584186476150168101766032692499585019572300009477107342325731806003539645749720998260074894005503243715777788187754550286488316966320623707403846973934258796072434720381188302956748572794547906591019552569722211941807674118070825610462581312695693650893380241875873254675082431417454526920871842555422162216221640850141209989864308"),
           e_float("0.081319696008915424923862092389736899811819511475112310887916890436246797704730704673729056088324784096921353721791055102563812592375189346816348952621603248827528923433998597558478366712658235656924601774921095079356565867874041683517103417362062033545480747089630436317004907748679995334911923789802669273406449877710501792905734235548100216242648355048174292232350538531258791292074224222773259366265"),
           e_float("0.083225259523836973901629476483826882570495097409722948211156651871411915926347698443317606324912095013846386144848076019334983060323996364575622231284156371311209803383868297666027142776492604713451401116304951916037351315288808149328656282592946622481248077426878390562084748274992526750980236994304445617883874613002677387890383751048690275480465558930319669630406765794141303531112767626885032234711"),
           e_float("0.085130518774003530041015433371012162123590338565297712899451572274044870858132401685294989726432756808251618198809497117031428754224182463768189333992610411750324249283786647605964593556580841488686189374310703328490638378071341189999989736919850819640255203388192531225871093915262442913255553964348789674386732423763999398384843950222967487918049677831350991752153610543052586794779289758903475438463"),
           e_float("0.087035466793942804442393380943589773874807384137701777807755858606555806090045313793349743667697666495784828540247587627381581545620918795107877715086227932192756144780086856515641298195558204484728115519445514732755572178929839231921546119143469647706999890259163924570037168105911942268212237224359201159897493303064534985539474450028748960081579369980283792804593122976296597796047722827878891322393"),
           e_float("0.088940096619320340510800454481897179643787714504988904779549430422905882044782054375828906100456879992313640144934242418129400679706865644316664458061953684856607341811667576418163339843372041295322709797007871604740491453013916791024918954258557372755707588171502529270934617054742841996068375532958173787119651672031238405892943090742861150287144368470297356054911834742041692715934237148148936307154"),
           e_float("0.090844401286964974994199780075024285173534240616972938308566124002318176156233343333797652771250228033597767340846476094210813299723066344955542455246897302945908336756464582189829335315055318509171623072382573656172991307091055717376996812421994002305994837615827430425975713541786272473174908117299253875064746541808768776888139494200307162476704515967971265236163183394620211890252072718610942571965"),
           e_float("0.092748373834894294768837248013614869545920040757666771889792264347367335922283987944081370356939214799510822558638935486112190359701790662792925773186897671244831971450474952036618937833706899006048469791213129922302798462311426294239282200248625552430083292818594731518350634587677158329361662244008606299358410273467156168469937529868888478684257396004386633746334360279683178784117372066095688934486"),
           e_float("0.094652007302340089278624856973167138217258137565762699413364163993950641989798336903224059425798872308451866930780842878435503832366316245159154044350517727312654980455290163514021860802028251981954724580292247887851185639950862672925642859475325430465650514258094310141033116303195908372018392349011340787373492402788283020594141795669947097771918542864972570255253465754192708165241608041715996375567"),
           e_float("0.09655529472977379853549858833033074225823562054271495313739665642376685099661084023094270272485976247900824483810911634635819558334630910267353320029261330296977292720266655308513559530586843550229208517388789783011887450865488554143475302590353915732321663418057567573042594801866258948380684000769091353165879953111046260532796891917772727185993569684246844052518121013717183610828519193371796413317"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00112_cos_x_near_pi_half(const bool b_write_output)
    {
      return TestCase_case_00112_cos_x_near_pi_half().execute(b_write_output);
    }
  }
}