data = readtable('data_upd3xlsx.xlsx');
POSCUMr = data.POSCUM;
NEGCUMr = data.NEGCUM;
DCUMr = data.DCUM;
Hr = data.H;

%replace all missing values (recorded as NaN's)  with 0's
POSCUMr(isnan(POSCUMr))=0;
NEGCUMr(isnan(NEGCUMr))=0;
DCUMr(isnan(DCUMr))=0;
Hr(isnan(Hr)) = 0;

%NOTE: Shift each variable one day forward, because the data is recorded at the end of each day
POSNr = [0;diff(POSCUMr(1:end-1))]; %# of NEW positive tests at t (raw)
NEGNr = [0;diff(NEGCUMr(1:end-1))]; %# of NEW positive tests at t (raw)
DNr = [0;diff(DCUMr(1:end-1))]; %# of NEW fatal cases at t (raw)

po = POSNr./(328.2*1e6);
ne = NEGNr./(328.2*1e6);
dn= DNr./(328.2*1e6);
h = Hr./(328.2*1e6);

daa = [POSNr NEGNr DNr Hr(1:end-1)];

doo = [po ne dn h(1:end-1)];