function [t,R2,RMSE] = tight_axs(ax,cb,x0,w0,mg,sims,varname,letter,labels,ref)
    num = {'~3,000,000','750,000','333,333','187,500','83,333','49,682','11,718','2,755'};
    %num = {'664,724','427,478','254,953','148,514','109,358','58,749','14,536','~3,000'};
    for i = 1 : 4
        % second row
        ax(1+i).Position(3) = w0;
        ax(1+i).Position(1) = x0 + (i-1)*(w0+mg);
        ax(1+i).Position(2) = ax(1+i).Position(2) + 0.04;
        % third row
        ax(5+i).Position(3) = w0;
        ax(5+i).Position(1) = x0 + (i-1)*(w0+mg);
        ax(5+i).Position(2) = ax(1+i).Position(2) - mg - ax(5+i).Position(4);
    end
    for i = 1 : 8
        add_title(ax(i+1),[letter{i} labels{i}],15,'in');
        if nargin == 9
            if i > 1
                [R2(i),RMSE(i),~,PBIAS] = estimate_evaluation_metric(sims(1).(varname),sims(i).(varname));
                R2(i)    = round(R2(i),2);
                RMSE(i)  = round(RMSE(i),2);
                str   = {['R^{2} = ' num2str(R2(i))], ['RMSE= ' num2str(RMSE(i))]};
                t(i)  = add_annot(ax(i+1),str,12,'lr');
            end
        else
            [R2(i),RMSE(i),~,PBIAS] = estimate_evaluation_metric(ref(1).(varname),sims(i).(varname));
            R2(i)    = round(R2(i),2);
            RMSE(i)  = round(RMSE(i),2);
            str   = {['R^{2} = ' num2str(R2(i))], ['RMSE= ' num2str(RMSE(i))]};
            t(i)  = add_annot(ax(i+1),str,12,'lr');
        end
        tt(i) = add_annot(ax(i+1),num{i},12,'ur');
        tt(i).Color = 'r';
    end

    cb.Position(1) = ax(6).Position(1);
    cb.Position(2) = ax(6).Position(2) - 0.075;
    cb.Position(3) = ax(9).Position(1) + ax(9).Position(3) - ax(6).Position(1);
    cb.Position(4) = 0.02;
    cb.FontSize = 14;
end