function [ E,conf_true ] = real_gml( fname )
%REAL_GML real gmal file and compute A and true_conf
%   Detailed explanation goes here
%Extracting edges from gml file graph
    inputfile = fopen(fname);
    A=[];
    l=0;
    k=1;
    id2num=containers.Map;
    id2color=containers.Map;
	value2num=containers.Map;
    edge_ids={};
    Ei=[];
    Ej=[];
    n=0;
    q=0;
    while 1
          tline = fgetl(inputfile); % Get a line from the input file
          if ~ischar(tline)% Quit if end of file
              break
          end

          %read node
          node = regexp(tline,'node','match');
          if ~isempty(node)
              while 1
                  line = fgetl(inputfile); % Get a line from the input file
                  assert( ischar(tline) );
                  ids=regexp( line,'id.+\d+','match');
                  lbs=regexp(line,'label.+','match');
                  if(~isempty(ids))
                      if(~isempty(lbs))
                          continue;
                      end
                      ids=regexp(line,'\d+','match');
                      id=ids{1};
                      %line
                      assert( ~isKey(id2num, id)); % id should not exist.
                      n = n+1;
                      id2num(id)=n; %map id to number
                  end
                  values = regexp( line,'value.+\d+','match');
                  if(~isempty(values))
                      values = regexp( line,'\d+','match');
                      value=values{1};
                      if( ~isKey(value2num,value))
                          q = q+1;
                          value2num(value)=q;
                      end
                      id2color(id)=value2num(value);
                  end
                  nend=regexp(line,']','match');
                  if(~isempty(nend))
                      break;
                  end
              end
          end
          
          %read edge
          edge = regexp(tline,'edge','match');
          if ~isempty(edge)
              while 1
                  line = fgetl(inputfile); % Get a line from the input file
                  assert( ischar(tline) );
                  sources=regexp( line,'source.+\d+','match');
                  if(~isempty(sources))
                      recordflag=true;
                      sources=regexp(line,'\d+','match');
                      source=sources{1};
                      if ( ~isKey(id2num, source) )
                          recordflag=false;
                          continue;
                      end
                      %assert( isKey(id2num, source) ); % id should exist.
                  end
                  targets=regexp( line,'target.+\d+','match');
                  if(~isempty(targets))
                      targets=regexp(line,'\d+','match');
                      target=targets{1};
                      %fprintf('%s -> %s\n',source, target);
                      if ( ~isKey(id2num, target) )
                          recordflag=false;
                          continue;
                      end
                      %assert( isKey(id2num, target) ); % id should exist.
                  end
                  nend=regexp(line,']','match');
                  if(~isempty(nend))
                      if(recordflag)
                          Ei=[Ei;id2num(source)];
                          Ej=[Ej;id2num(target)];
                      end
                      break;
                  end
              end
          end
    end
    E=[Ei Ej];
    conf_true=zeros(n,1);
    %id2color.keys()
    for i=id2color.keys()
        a=id2num(i{1});
        q=id2color(i{1});
        %[a q]
        conf_true(a,1)=q;
    end
    %A=E2A(E);
    %E=A2E(A);
    %id2num
    %conf_true

end

