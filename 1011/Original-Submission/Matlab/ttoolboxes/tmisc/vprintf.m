function varargout = vprintf(varargin); 
% [ ret ] = vprintf( string, args, [options])
% vprintf(string, args, [options], 'save')
% [ ret ] = vprintf( string, 'load', ['str', str] )
%
% More powerful version of sprintf.
%
% Input: 
% %%%%%%%%%%%
%   string              A format string as used by fprintf
%                       vprintf knows the special format strings %v and %r
%                           '%v'        passes argument to vdisp, i.e. prints out nearly anything
%                           '%r'        very compact printing of arrays
%   args                The arguments as given for fprintf
%
%
% Options: 
% %%%%%%%%%%%%%%
%   'cpr',val           string, Uses cprintf to print out the text. Uses the format defined by val
%   'imp',[imp val]     Prints only if imp>=val. If imp>=val-1 and 'string',str is given, then the output is appended to str.
%                           if val>=importance   : varargin is passed to fprintf/cprintf
%                           if val>=importance-1 : if str is given (as a string), then the output is appended to str and returned
%                           if val<importance-1  : nothing happens
%   'sze',[imp val]     Similar behaviour as for imp. But does print 'BIG SIZE SKIP\n' instead of not printing
%   'str',str           string, default='', Output is appended to str and returned
%   'npr'               No print output is generated
%   'once',val          Only executes the command once. Possible values:
%                           val==1  : saves the command such that further calls of it are not exectued
%                           val==-1 : Clears all saved commands, i.e. all commands are executed again when called
%                           E.g.: vprintf('test\n','once',1); vprintf('test\n','once',1); vprintf('once',-1); vprintf('test\n','once',1);
%                       
%   'save'              (experimental) Does not print the message, until vprintf('load') is called.
%   'load'              (experimental) Prints out all saved messages
%
%
% Output
% %%%%%%%%%%%
%   ret                 str, the printed string
%                       Not set if the options 'save' is used
%
% Note:
% %%%%%%%
%  Options'save' and 'once' store all input in a persistent variable. Thus, do not use these options if you pass much data.
%
%
% E.g.: [str]=vprintf('Talk %i me: \n%v\n',2,[1 2 3; 3 2 1],'cpr',[.5 .5 0],'str','')
%
% See also: vdisp, fprintf, cprintf
%
% Written by: tommsch, 2018

%Changelog: 2019_10_08 - tommsch - Added option 'npr' %XX Still needs to be tested
%           2019_01_08 - tommsch - Added option 'once'
%           2020_02_05 - tommsch - Added options 'save' and 'load'

% Add option such that output is only in return value/screen/nowhere/both
% Maybe remove behaviour that output in 'ret' is different to screen output (except for option 'npr')

%#ok<*ALIGN>



persistent oncearray;
persistent savemessage;

[saveflag,varargin] = parsem( {'save'},varargin );
if(saveflag)
    savemessage{end+1} = varargin; 
    return; end;
[loadflag,varargin] = parsem( {'load'}, varargin );
if(loadflag)
    [str,varargin] = parsem( {'str'}, varargin, '' );
    parsem( 'test', varargin );
    for i = 1:numel( savemessage )
        str = vprintf( savemessage{i}{:}, 'str',str ); end;
    savemessage = {}; 
    varargout{1} = str;
    return; end
[oncearray_flag,varargin] = parsem( {'once'}, varargin, 0 );
if( isempty(oncearray_flag) )
    oncearray = {}; end; 
if( oncearray_flag==-1 ); %reset
    oncearray = {}; 
    return; end;
if( oncearray_flag==1 )
    if( searchincellarray(varargin,oncearray,0) )
        return; 
    else
        oncearray{end+1} = varargin; end; end;




[str,varargin] = parsem( {'str'}, varargin, '' );
[imp,varargin] = parsem( {'imp'}, varargin, [1 1] );
if( imp(2)<imp(1)-1 ); 
    if( nargout==1 ); 
        varargout{1} = str; end;
    return; end;
if( imp(2)<imp(1) && nargout==0 ); 
    if( nargout==1 ); 
        varargout{1} = str; end;
    return; end;
[npr,varargin] = parsem( {'noprint','npr'}, varargin );
[sze,varargin] = parsem( {'sze'}, varargin, [1 1] );
if( sze(2)>sze(1)*4-1 ); 
    txt = sprintf('BIG SIZE SKIP.');
    if( imp(2)>=imp(1) && ~npr ); 
        disp(txt); end;
    if( nargout==1 ); 
        varargout{1} = [str txt]; end;
    return; end;
[format,varargin] = parsem( {'cpr'}, varargin, '' );



%parse input string and replace '%v' by '%s%' and the corresponding argument by the output of vdisp(argument)
percentidx = strfind( varargin{1}, '%' );
i = 1; %counts the pairs "%-something"
j = 1; %counts the numbers in percentidx
while( true )
    if( j>size(percentidx,2) ); 
        break; end;
    if( varargin{1}(percentidx(j)+1)=='%' ); 
        j=j+2;
    elseif( varargin{1}(percentidx(j)+1)~='v' && varargin{1}(percentidx(j)+1)~='r' ); 
        i=i+1; j=j+1;
    elseif( varargin{1}(percentidx(j)+1)=='v' ); 
        varargin{i+1} = vdisp(varargin{1+i});
        varargin{1} = [varargin{1}(1:percentidx(j)-1) '%s' varargin{1}(percentidx(j)+2:end)];
        i = i+1; j = j+1; 
    elseif( varargin{1}(percentidx(j)+1)=='r' ); 
        varargin{i+1} = sprintf( '%i ', varargin{i+1} );
        varargin{1} = [varargin{1}(1:percentidx(j)-1) '%s' varargin{1}(percentidx(j)+2:end)];
        i = i+1; 
        j = j+1; end; end;

if( size(varargin,2)~=i );  %test if number of arguments is correct
    error( 'vprintf:nargin', 'Number of arguments wrong' ); end;


if( imp(2)>=imp(1) && ~npr )
    if( isempty(format) )
        fprintf( varargin{:} );  
    else
        cprintf( format, varargin{:} ); end; end

if( nargout==1 ); 
    txt = sprintf( varargin{:} );
    varargout{1} = [str txt]; end;

end

function count = cprintf( style, format, varargin )
% CPRINTF displays styled formatted text in the Command Window
%
% Syntax:
%    count = cprintf(style,format,...)
%
% Description:
%    CPRINTF processes the specified text using the exact same FORMAT
%    arguments accepted by the built-in SPRINTF and FPRINTF functions.
%
%    CPRINTF then displays the text in the Command Window using the
%    specified STYLE argument. The accepted styles are those used for
%    Matlab's syntax highlighting (see: File / Preferences / Colors / 
%    M-file Syntax Highlighting Colors), and also user-defined colors.
%
%    The possible pre-defined STYLE names are:
%
%       'Text'                 - default: black
%       'Keywords'             - default: blue
%       'Comments'             - default: green
%       'Strings'              - default: purple
%       'UnterminatedStrings'  - default: dark red
%       'SystemCommands'       - default: orange
%       'Errors'               - default: light red
%       'Hyperlinks'           - default: underlined blue
%
%       'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White'
%
%    STYLE beginning with '-' or '_' will be underlined. For example:
%          '-Blue' is underlined blue, like 'Hyperlinks';
%          '_Comments' is underlined green etc.
%
%    STYLE beginning with '*' will be bold (R2011b+ only). For example:
%          '*Blue' is bold blue;
%          '*Comments' is bold green etc.
%    Note: Matlab does not currently support both bold and underline,
%          only one of them can be used in a single cprintf command. But of
%          course bold and underline can be mixed by using separate commands.
%
%    STYLE also accepts a regular Matlab RGB vector, that can be underlined
%    and bolded: -[0,1,1] means underlined cyan, '*[1,0,0]' is bold red.
%
%    STYLE is case-insensitive and accepts unique partial strings just
%    like handle property names.
%
%    CPRINTF by itself, without any input parameters, displays a demo
%
% Example:
%    cprintf;   % displays the demo
%    cprintf('text',   'regular black text');
%    cprintf('hyper',  'followed %s','by');
%    cprintf('key',    '%d colored', 4);
%    cprintf('-comment','& underlined');
%    cprintf('err',    'elements\n');
%    cprintf('cyan',   'cyan');
%    cprintf('_green', 'underlined green');
%    cprintf(-[1,0,1], 'underlined magenta');
%    cprintf([1,0.5,0],'and multi-\nline orange\n');
%    cprintf('*blue',  'and *bold* (R2011b+ only)\n');
%    cprintf('string');  % same as fprintf('string') and cprintf('text','string')
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab
%    functionality. It works on Matlab 7+, but use at your own risk!
%
%    A technical description of the implementation can be found at:
%    <a href="http://undocumentedmatlab.com/blog/cprintf/">http://UndocumentedMatlab.com/blog/cprintf/</a>
%
% Limitations:
%    1. In R2011a and earlier, a single space char is inserted at the
%       beginning of each CPRINTF text segment (this is ok in R2011b+).
%
%    2. In R2011a and earlier, consecutive differently-colored multi-line
%       CPRINTFs sometimes display incorrectly on the bottom line.
%       As far as I could tell this is due to a Matlab bug. Examples:
%         >> cprintf('-str','under\nline'); cprintf('err','red\n'); % hidden 'red', unhidden '_'
%         >> cprintf('str','regu\nlar'); cprintf('err','red\n'); % underline red (not purple) 'lar'
%
%    3. Sometimes, non newline ('\n')-terminated segments display unstyled
%       (black) when the command prompt chevron ('>>') regains focus on the
%       continuation of that line (I can't pinpoint when this happens). 
%       To fix this, simply newline-terminate all command-prompt messages.
%
%    4. In R2011b and later, the above errors appear to be fixed. However,
%       the last character of an underlined segment is not underlined for
%       some unknown reason (add an extra space character to make it look better)
%
%    5. In old Matlab versions (e.g., Matlab 7.1 R14), multi-line styles
%       only affect the first line. Single-line styles work as expected.
%       R14 also appends a single space after underlined segments.
%
%    6. Bold style is only supported on R2011b+, and cannot also be underlined.
%
% Change log:
%    2015-06-24: Fixed a few discoloration issues (some other issues still remain)
%    2015-03-20: Fix: if command window isn't defined yet (startup) use standard fprintf as suggested by John Marozas
%    2012-08-09: Graceful degradation support for deployed (compiled) and non-desktop applications; minor bug fixes
%    2012-08-06: Fixes for R2012b; added bold style; accept RGB string (non-numeric) style
%    2011-11-27: Fixes for R2011b
%    2011-08-29: Fix by Danilo (FEX comment) for non-default text colors
%    2011-03-04: Performance improvement
%    2010-06-27: Fix for R2010a/b; fixed edge case reported by Sharron; CPRINTF with no args runs the demo
%    2009-09-28: Fixed edge-case problem reported by Swagat K
%    2009-05-28: corrected nargout behavior suggested by Andreas Gäb
%    2009-05-13: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
%
% See also:
%    sprintf, fprintf
%
% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.
%
% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.10 $  $Date: 2015/06/24 01:29:18 $

  persistent majorVersion minorVersion
  if isempty(majorVersion)
      %v = version; if str2double(v(1:3)) <= 7.1
      %majorVersion = str2double(regexprep(version,'^(\d+).*','$1'));
      %minorVersion = str2double(regexprep(version,'^\d+\.(\d+).*','$1'));
      %[a,b,c,d,versionIdStrs]=regexp(version,'^(\d+)\.(\d+).*');  %#ok unused
      v = sscanf(version, '%d.', 2);
      majorVersion = v(1); %str2double(versionIdStrs{1}{1});
      minorVersion = v(2); %str2double(versionIdStrs{1}{2});
  end

  % The following is for debug use only:
  %global docElement txt el
  if ~exist('el','var') || isempty(el),  el=handle([]);  end  %#ok mlint short-circuit error ("used before defined")
  %if nargin<1, showDemo(majorVersion,minorVersion); return;  end
  if isempty(style),  return;  end
  if all(ishandle(style)) && length(style)~=3
      dumpElement(style);
      return;
  end

  % Process the text string
  if nargin<2, format = style; style='text';  end
  %error(nargchk(2, inf, nargin, 'struct'));
  %str = sprintf(format,varargin{:});

  % In compiled mode
  try useDesktop = usejava('desktop'); catch, useDesktop = false; end
  if isdeployed | ~useDesktop %#ok<OR2> - for Matlab 6 compatibility
      % do not display any formatting - use simple fprintf()
      % See: http://undocumentedmatlab.com/blog/bold-color-text-in-the-command-window/#comment-103035
      % Also see: https://mail.google.com/mail/u/0/?ui=2&shva=1#all/1390a26e7ef4aa4d
      % Also see: https://mail.google.com/mail/u/0/?ui=2&shva=1#all/13a6ed3223333b21
      count1 = fprintf(format,varargin{:});
  else
      % Else (Matlab desktop mode)
      % Get the normalized style name and underlining flag
      [underlineFlag, boldFlag, style, debugFlag] = processStyleInfo(style);

      % Set hyperlinking, if so requested
      if underlineFlag
          format = ['<a href="">' format '</a>'];

          % Matlab 7.1 R14 (possibly a few newer versions as well?)
          % have a bug in rendering consecutive hyperlinks
          % This is fixed by appending a single non-linked space
          if majorVersion < 7 || (majorVersion==7 && minorVersion <= 1)
              format(end+1) = ' ';
          end
      end

      % Set bold, if requested and supported (R2011b+)
      if boldFlag
          if (majorVersion > 7 || minorVersion >= 13)
              format = ['<strong>' format '</strong>'];
          else
              boldFlag = 0;
          end
      end

      % Get the current CW position
      cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
      lastPos = cmdWinDoc.getLength;

      % If not beginning of line
      bolFlag = 0;  %#ok
      %if docElement.getEndOffset - docElement.getStartOffset > 1
          % Display a hyperlink element in order to force element separation
          % (otherwise adjacent elements on the same line will be merged)
          if majorVersion<7 || (majorVersion==7 && minorVersion<13)
              if ~underlineFlag
                  fprintf('<a href=""> </a>');  %fprintf('<a href=""> </a>\b');
              elseif format(end)~=10  % if no newline at end
                  fprintf(' ');  %fprintf(' \b');
              end
          end
          %drawnow;
          bolFlag = 1;
      %end

      % Get a handle to the Command Window component
      mde = com.mathworks.mde.desk.MLDesktop.getInstance;
      cw = mde.getClient('Command Window');

      % Fix: if command window isn't defined yet (startup), use standard fprintf()
      if (isempty(cw))
         count1 = fprintf(format,varargin{:});
         if nargout
             count = count1;
         end
         return;
      end
      
      xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);

      % Store the CW background color as a special color pref
      % This way, if the CW bg color changes (via File/Preferences), 
      % it will also affect existing rendered strs
      com.mathworks.services.Prefs.setColorPref('CW_BG_Color',xCmdWndView.getBackground);

      % Display the text in the Command Window
      % Note: fprintf(2,...) is required in order to add formatting tokens, which
      % ^^^^  can then be updated below (no such tokens when outputting to stdout)
      count1 = fprintf(2,format,varargin{:});

      % Repaint the command window
      %awtinvoke(cmdWinDoc,'remove',lastPos,1);   % TODO: find out how to remove the extra '_'
      drawnow;  % this is necessary for the following to work properly (refer to Evgeny Pr in FEX comment 16/1/2011)
      xCmdWndView.repaint;
      %hListeners = cmdWinDoc.getDocumentListeners; for idx=1:numel(hListeners), try hListeners(idx).repaint; catch, end, end

      docElement = cmdWinDoc.getParagraphElement(lastPos+1);
      if majorVersion<7 || (majorVersion==7 && minorVersion<13)
          if bolFlag && ~underlineFlag
              % Set the leading hyperlink space character ('_') to the bg color, effectively hiding it
              % Note: old Matlab versions have a bug in hyperlinks that need to be accounted for...
              %disp(' '); dumpElement(docElement)
              setElementStyle(docElement,'CW_BG_Color',1+underlineFlag,majorVersion,minorVersion); %+getUrlsFix(docElement));
              %disp(' '); dumpElement(docElement)
              el(end+1) = handle(docElement);  % #ok used in debug only
          end

          % Fix a problem with some hidden hyperlinks becoming unhidden...
          fixHyperlink(docElement);
          %dumpElement(docElement);
      end

      % Get the Document Element(s) corresponding to the latest fprintf operation
      while docElement.getStartOffset < cmdWinDoc.getLength
          % Set the element style according to the current style
          if debugFlag, dumpElement(docElement); end
          specialFlag = underlineFlag | boldFlag;
          setElementStyle(docElement,style,specialFlag,majorVersion,minorVersion);
          if debugFlag, dumpElement(docElement); end
          docElement2 = cmdWinDoc.getParagraphElement(docElement.getEndOffset+1);
          if isequal(docElement,docElement2),  break;  end
          docElement = docElement2;
      end
      if debugFlag, dumpElement(docElement); end

      % Force a Command-Window repaint
      % Note: this is important in case the rendered str was not '\n'-terminated
      xCmdWndView.repaint;

      % The following is for debug use only:
      el(end+1) = handle(docElement);  %#ok used in debug only
      %elementStart  = docElement.getStartOffset;
      %elementLength = docElement.getEndOffset - elementStart;
      %txt = cmdWinDoc.getText(elementStart,elementLength);
  end

  if nargout
      count = count1;
  end
  return;  % debug breakpoint
  
end

% Process the requested style information
function [ underlineFlag, boldFlag, style, debugFlag ] = processStyleInfo( style )
  underlineFlag = 0;
  boldFlag = 0;
  debugFlag = 0;

  % First, strip out the underline/bold markers
  if( ischar(style) )
      % Styles containing '-' or '_' should be underlined (using a no-target hyperlink hack)
      %if style(1)=='-'
      underlineIdx = (style=='-') | (style=='_');
      if( any(underlineIdx) )
          underlineFlag = 1;
          %style = style(2:end);
          style = style(~underlineIdx); end;

      % Check for bold style (only if not underlined)
      boldIdx = (style=='*');
      if( any(boldIdx) )
          boldFlag = 1;
          style = style(~boldIdx); end;
      if( underlineFlag && boldFlag )
          warning('YMA:cprintf:BoldUnderline','Matlab does not support both bold & underline'); end;

      % Check for debug mode (style contains '!')
      debugIdx = (style=='!');
      if( any(debugIdx) )
          debugFlag = 1;
          style = style(~debugIdx); end;

      % Check if the remaining style sting is a numeric vector
      %styleNum = str2num(style); %#ok<ST2NM>  % not good because style='text' is evaled!
      %if ~isempty(styleNum)
      if( any(style==' ' | style==',' | style==';') )
          style = str2num(style); end; %#ok<ST2NM>
  end

  % Style = valid matlab RGB vector
  if isnumeric(style) && length(style)==3 && all(style<=1) && all(abs(style)>=0)
      if any(style<0)
          underlineFlag = 1;
          style = abs(style);
      end
      style = getColorStyle(style);

  elseif ~ischar(style)
      error('YMA:cprintf:InvalidStyle','Invalid style - see help section for a list of valid style values')

  % Style name
  else
      % Try case-insensitive partial/full match with the accepted style names
      matlabStyles = {'Text','Keywords','Comments','Strings','UnterminatedStrings','SystemCommands','Errors'};
      validStyles  = [matlabStyles, ...
                      'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White', ...
                      'Hyperlinks'];
      matches = find(strncmpi(style,validStyles,length(style)));

      % No match - error
      if isempty(matches)
          error('YMA:cprintf:InvalidStyle','Invalid style - see help section for a list of valid style values')

      % Too many matches (ambiguous) - error
      elseif length(matches) > 1
          error('YMA:cprintf:AmbigStyle','Ambiguous style name - supply extra characters for uniqueness')

      % Regular text
      elseif matches == 1
          style = 'ColorsText';  % fixed by Danilo, 29/8/2011

      % Highlight preference style name
      elseif matches <= length(matlabStyles)
          style = ['Colors_M_' validStyles{matches}];

      % Color name
      elseif matches < length(validStyles)
          colors = [0,0,0; 0,1,1; 1,0,1; 0,0,1; 0,1,0; 1,0,0; 1,1,0; 1,1,1];
          requestedColor = colors(matches-length(matlabStyles),:);
          style = getColorStyle(requestedColor);

      % Hyperlink
      else
          style = 'Colors_HTML_HTMLLinks';  % CWLink
          underlineFlag = 1;
      end
  end
end

% Convert a Matlab RGB vector into a known style name (e.g., '[255,37,0]')
function styleName = getColorStyle( rgb )
  intColor = int32( rgb*255 );
  javaColor = java.awt.Color( intColor(1), intColor(2), intColor(3) );
  styleName = sprintf( '[%d,%d,%d]', intColor );
  com.mathworks.services.Prefs.setColorPref( styleName, javaColor );
end

% Fix a bug in some Matlab versions, where the number of URL segments
% is larger than the number of style segments in a doc element
function delta = getUrlsFix(docElement)  %#ok currently unused
  tokens = docElement.getAttribute('SyntaxTokens');
  links  = docElement.getAttribute('LinkStartTokens');
  if length(links) > length(tokens(1))
      delta = length(links) > length(tokens(1));
  else
      delta = 0;
  end
end

% fprintf(2,str) causes all previous '_'s in the line to become red - fix this
function fixHyperlink(docElement)
  try
      tokens = docElement.getAttribute('SyntaxTokens');
      urls   = docElement.getAttribute('HtmlLink');
      urls   = urls(2);
      links  = docElement.getAttribute('LinkStartTokens');
      offsets = tokens(1);
      styles  = tokens(2);
      doc = docElement.getDocument;

      % Loop over all segments in this docElement
      for idx = 1 : length(offsets)-1
          % If this is a hyperlink with no URL target and starts with ' ' and is collored as an error (red)...
          if strcmp(styles(idx).char,'Colors_M_Errors')
              character = char(doc.getText(offsets(idx)+docElement.getStartOffset,1));
              if strcmp(character,' ')
                  if isempty(urls(idx)) && links(idx)==0
                      % Revert the style color to the CW background color (i.e., hide it!)
                      styles(idx) = java.lang.String('CW_BG_Color');
                  end
              end
          end
      end
  catch
      % never mind...
  end
end

% Set an element to a particular style (color)
function setElementStyle(docElement,style,specialFlag, majorVersion,minorVersion)
  %global tokens links urls urlTargets  % for debug only
  global oldStyles
  if nargin<3,  specialFlag=0;  end
  % Set the last Element token to the requested style:
  % Colors:
  tokens = docElement.getAttribute('SyntaxTokens');
  try
      styles = tokens(2);
      oldStyles{end+1} = cell(styles);

      % Correct edge case problem
      extraInd = double(majorVersion>7 || (majorVersion==7 && minorVersion>=13));  % =0 for R2011a-, =1 for R2011b+
      %{
      if ~strcmp('CWLink',char(styles(end-hyperlinkFlag))) && ...
          strcmp('CWLink',char(styles(end-hyperlinkFlag-1)))
         extraInd = 0;%1;
      end
      hyperlinkFlag = ~isempty(strmatch('CWLink',tokens(2)));
      hyperlinkFlag = 0 + any(cellfun(@(c)(~isempty(c)&&strcmp(c,'CWLink')),cell(tokens(2))));
      %}

      jStyle = java.lang.String(style);
      if numel(styles)==4 && isempty(char(styles(2)))
          % Attempt to fix discoloration issues - NOT SURE THAT THIS IS OK! - 24/6/2015
          styles(1) = jStyle;
      end
      styles(end-extraInd) = java.lang.String('');
      styles(end-extraInd-specialFlag) = jStyle;  % #ok apparently unused but in reality used by Java
      if extraInd
          styles(end-specialFlag) = jStyle;
      end

      oldStyles{end} = [oldStyles{end} cell(styles)];
  catch
      % never mind for now
  end
  
  % Underlines (hyperlinks):
  %{
  links = docElement.getAttribute('LinkStartTokens');
  if isempty(links)
      %docElement.addAttribute('LinkStartTokens',repmat(int32(-1),length(tokens(2)),1));
  else
      %TODO: remove hyperlink by setting the value to -1
  end
  %}

  % Correct empty URLs to be un-hyperlinkable (only underlined)
  urls = docElement.getAttribute('HtmlLink');
  if ~isempty(urls)
      urlTargets = urls(2);
      for urlIdx = 1 : length(urlTargets)
          try
              if urlTargets(urlIdx).length < 1
                  urlTargets(urlIdx) = [];  % '' => []
              end
          catch
              % never mind...
              a=1;  %#ok used for debug breakpoint...
          end
      end
  end
  
  % Bold: (currently unused because we cannot modify this immutable int32 numeric array)
  %{
  try
      %hasBold = docElement.isDefined('BoldStartTokens');
      bolds = docElement.getAttribute('BoldStartTokens');
      if ~isempty(bolds)
          %docElement.addAttribute('BoldStartTokens',repmat(int32(1),length(bolds),1));
      end
  catch
      % never mind - ignore...
      a=1;  %#ok used for debug breakpoint...
  end
  %}
  
  return;  % debug breakpoint
end

% Display information about element(s)
function dumpElement(docElements)
  %return;
  disp(' ');
  numElements = length(docElements);
  cmdWinDoc = docElements(1).getDocument;
  for elementIdx = 1 : numElements
      if numElements > 1,  fprintf('Element #%d:\n',elementIdx);  end
      docElement = docElements(elementIdx);
      if ~isjava(docElement),  docElement = docElement.java;  end
      %docElement.dump(java.lang.System.out,1)
      disp(docElement)
      tokens = docElement.getAttribute('SyntaxTokens');
      if isempty(tokens),  continue;  end
      links = docElement.getAttribute('LinkStartTokens');
      urls  = docElement.getAttribute('HtmlLink');
      try bolds = docElement.getAttribute('BoldStartTokens'); catch, bolds = []; end
      txt = {};
      tokenLengths = tokens(1);
      for tokenIdx = 1 : length(tokenLengths)-1
          tokenLength = diff(tokenLengths(tokenIdx+[0,1]));
          if (tokenLength < 0)
              tokenLength = docElement.getEndOffset - docElement.getStartOffset - tokenLengths(tokenIdx);
          end
          txt{tokenIdx} = cmdWinDoc.getText(docElement.getStartOffset+tokenLengths(tokenIdx),tokenLength).char;  %#ok
      end
      lastTokenStartOffset = docElement.getStartOffset + tokenLengths(end);
      try
          txt{end+1} = cmdWinDoc.getText(lastTokenStartOffset, docElement.getEndOffset-lastTokenStartOffset).char; %#ok
      catch
          txt{end+1} = ''; %#ok<AGROW>
      end
      %cmdWinDoc.uiinspect
      %docElement.uiinspect
      txt = strrep(txt',sprintf('\n'),'\n'); %#ok<SPRINTFN>
      try
          data = [cell(tokens(2)) m2c(tokens(1)) m2c(links) m2c(urls(1)) cell(urls(2)) m2c(bolds) txt];
          if elementIdx==1
              disp('    SyntaxTokens(2,1) - LinkStartTokens - HtmlLink(1,2) - BoldStartTokens - txt');
              disp('    ==============================================================================');
          end
      catch
          try
              data = [cell(tokens(2)) m2c(tokens(1)) m2c(links) txt];
          catch
              disp([cell(tokens(2)) m2c(tokens(1)) txt]);
              try
                  data = [m2c(links) m2c(urls(1)) cell(urls(2))];
              catch
                  % Mtlab 7.1 only has urls(1)...
                  data = [m2c(links) cell(urls)];
              end
          end
      end
      disp(data)
  end
end

% Utility function to convert matrix => cell
function cells = m2c( data )
  %datasize = size(data);  cells = mat2cell(data,ones(1,datasize(1)),ones(1,datasize(2)));
  cells = num2cell( data );
end
% 
% % Display the help and demo
% function showDemo(majorVersion,minorVersion)
%   fprintf('cprintf displays formatted text in the Command Window.\n\n');
%   fprintf('Syntax: count = cprintf(style,format,...);  click <a href="matlab:help cprintf">here</a> for details.\n\n');
%   url = 'http://UndocumentedMatlab.com/blog/cprintf/';
%   fprintf(['Technical description: <a href="' url '">' url '</a>\n\n']);
%   fprintf('Demo:\n\n');
%   boldFlag = majorVersion>7 || (majorVersion==7 && minorVersion>=13);
%   s = ['cprintf(''text'',    ''regular black text'');' 10 ...
%        'cprintf(''hyper'',   ''followed %s'',''by'');' 10 ...
%        'cprintf(''key'',     ''%d colored'',' num2str(4+boldFlag) ');' 10 ...
%        'cprintf(''-comment'',''& underlined'');' 10 ...
%        'cprintf(''err'',     ''elements:\n'');' 10 ...
%        'cprintf(''cyan'',    ''cyan'');' 10 ...
%        'cprintf(''_green'',  ''underlined green'');' 10 ...
%        'cprintf(-[1,0,1],  ''underlined magenta'');' 10 ...
%        'cprintf([1,0.5,0], ''and multi-\nline orange\n'');' 10];
%    if boldFlag
%        % In R2011b+ the internal bug that causes the need for an extra space
%        % is apparently fixed, so we must insert the sparator spaces manually...
%        % On the other hand, 2011b enables *bold* format
%        s = [s 'cprintf(''*blue'',   ''and *bold* (R2011b+ only)\n'');' 10];
%        s = strrep(s, ''')',' '')');
%        s = strrep(s, ''',5)',' '',5)');
%        s = strrep(s, '\n ','\n');
%    end
%    disp(s);
%    eval(s);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%
% - Fix: Remove leading space char (hidden underline '_')
% - Fix: Find workaround for multi-line quirks/limitations
% - Fix: Non-\n-terminated segments are displayed as black
% - Fix: Check whether the hyperlink fix for 7.1 is also needed on 7.2 etc.
% - Enh: Add font support



function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   