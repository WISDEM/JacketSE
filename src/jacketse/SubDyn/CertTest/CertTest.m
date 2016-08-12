
lgndTxt={'NREL results','CertTest'};

% data for looking at summary file results:
data = cell(2,1);
w1=3;    %width of first line
w2=1.5;  %width of second line
ls='--'; %line style


for test =1 %1:5
    
    RootName = sprintf ('Test%02.0f', test); 
    FileName = [ RootName '.SD.out' ];
    Files={[RootName filesep 'NREL_Results' filesep FileName],[RootName filesep FileName]};
    
    PlotFASToutput( Files, lgndTxt );
                
%%
                
    % check the summary file data:
    for k=1:2
        data{k} = ReadSubDynSummary(strrep(Files{k},'.out','.sum'));
    end
    
%%
    figure;
    subplot(2,2,1); hold on; ylabel('Node coords')
    plot( data{1}.Node_coords(:), 'LineWidth',w1 ); plot( data{2}.Node_coords(:), ls,'LineWidth',w2  ); legend(lgndTxt);
    
    subplot(2,2,2); hold on; ylabel('Mesh node')
    plot( data{1}.MeshNode, 'LineWidth',w1 ); plot( data{2}.MeshNode, ls,'LineWidth',w2  ); legend(lgndTxt);
    
    subplot(2,2,3); hold on; ylabel('CM coords')
    plot( data{1}.CM_coords, 'LineWidth',w1 ); plot( data{2}.CM_coords, ls,'LineWidth',w2  ); legend(lgndTxt);
    
    subplot(2,2,4); hold on; ylabel('CM')
    bar( [data{1}.Mass; data{2}.Mass]  ); set(gca,'xTick',[1 2], 'xTickLabel',lgndTxt);

    
    % note that this will have severe consequences if the sizes of these
    % values are not the same between NREL Results and new results:   
    figure;
    subplot(3,1,1); hold on; ylabel('FEM eigenvalues')
    bar( [data{1}.FEM_eig data{2}.FEM_eig]  ); legend(lgndTxt);
    
    subplot(3,1,2); hold on; ylabel('CB eigenvalues')
    bar( [data{1}.CB_eig data{2}.CB_eig]  ); legend(lgndTxt);
    
    subplot(3,1,3); hold on; ylabel('FEM eigenvector entries')
    plot( data{1}.FEM_eigVec(:), 'LineWidth',w1 ); plot( data{2}.FEM_eigVec(:), ls,'LineWidth',w2  ); legend(lgndTxt);
    
    
    figure;
    subplot(2,1,1); hold on; ylabel('Phi M')
    plot( data{1}.PhiM(:), 'LineWidth',w1 ); plot( data{2}.PhiM(:), ls,'LineWidth',w2  ); legend(lgndTxt);    

    subplot(2,1,2); hold on; ylabel('Phi R')
    plot( data{1}.PhiR(:), 'LineWidth',w1 ); plot( data{2}.PhiR(:), ls,'LineWidth',w2  ); legend(lgndTxt);

    
    figure;
    subplot(3,1,1); hold on; ylabel('KBBt')
    plot( data{1}.KBBt(:), 'LineWidth',w1 ); plot( data{2}.KBBt(:), ls, 'LineWidth',w2  ); legend(lgndTxt);    

    subplot(3,1,2); hold on; ylabel('MBBt')
    plot( data{1}.MBBt(:), 'LineWidth',w1 ); plot( data{2}.MBBt(:), ls,'LineWidth',w2  ); legend(lgndTxt);

    subplot(3,1,3); hold on; ylabel('MRB')
    plot( data{1}.MRB(:), 'LineWidth',w1 ); plot( data{2}.MRB(:), ls,'LineWidth',w2  ); legend(lgndTxt);
                           
                
end