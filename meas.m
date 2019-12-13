function qmes = meas(q)

qdeg = q * 180 / pi;
qmes = round(qdeg * 10) / 10;
qmes = qmes * pi / 180;
%qmes = q;