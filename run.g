Read("read.g");
SetInfoLevel(InfoCaratCat, 1);
SetInfoLevel(InfoCIB, 1);
CIBRunJob := function(queue, set, key)
    local sgrp, name, max, len, exp, cib, force, cibs;

    force := ValueOption("force")=true;

    while true do
        name := CIBFetchName(queue, set);
        if name = fail then
            QuitGap(0);
        fi;
        if not force then
            cibs := RedisCommand(StringFormatted("HGET {} {}", name, key));
            if cibs <> "" then
                CIBRemoveName( name, set );
                continue;
            fi;
        fi;
        sgrp := CaratCatFetchAffineCrystGroup(name : calculate:=true); # : subgroups:="cib1" );
        cib := CIB.CofiniteIntegralBraceVectorSystemsByContext( sgrp );
        SetCofiniteIntegralBraceVectorSystems( sgrp, cib );
        CofiniteIntegralBracesRepsGenerators( sgrp );
        CaratCatSaveCIBCocycles(sgrp, key );
        CIBRemoveName( name, set );
        # RedisCommand(StringFormatted("SREM {} {}", set, name));
        Unbind(sgrp);
    od;
end;
CIBRunJob1 := function(queue, set, key)
    local sgrp, name, max, len, exp, cib, force, cibs;

    force := ValueOption("force")=true;

    while true do
        name := CIBFetchName(queue, set);
        if name = fail then
            QuitGap(0);
        fi;
        if not force then
            cibs := RedisCommand(StringFormatted("HGET {} {}", name, key));
            if cibs <> "" then
                CIBRemoveName( name, set );
                continue;
            fi;
        fi;
        sgrp := CaratCatCofiniteIntegralBraceVectorSystems( name, key );
        CofiniteIntegralBracesRepsGenerators( sgrp );
        CaratCatSaveCIBCocycles(sgrp, key );
        CIBRemoveName( name, set );
        #RedisCommand(StringFormatted("SREM {} {}", set, name));
        Unbind(sgrp);
    od;
end;
#CIBRunJob("input", "running", "cib:f" : force:=true);

CIBRunJobCheck := function(queue, set, max_iters)
    local name, context, err, cnt;

    context := ValueOption("context")=true;

    cnt := 0;
    while true do
        cnt := cnt+1;
        if cnt = max_iters then
            Info(InfoCaratCat, 1, "Maximum number of iterations reached, quitting ...");
            QuitGap(1);
        fi;
        name := CIBFetchName(queue, set);
        if name = fail then
            QuitGap(0);
        fi;
        err := CaratCatRecalculate(name : context:=context);
        if not err then
            RedisCommand("SADD error {}", name);
        fi;
        CIBRemoveName( name, set );
        Info( InfoCaratCat, 1, cnt, " - ", name, ": ", err);
    od;
end;
CIBRunJobCheck( "ainput", "arunning", 3 );
