function [q_hat, bias_hat, P_hat] = mekf(...
		current_gyro, current_star_tracker_dcm, q_prev, bias_prev, P_prev, dt, params)
	% Prediction: undo the known gyro mounting/calibration transformation.
	%girdiler
	%Güncel gyro ölçümü
	%Startrackerden gelen DCM
	%Bir önceki quaternion
	%Önceki gyro bias kestirimi
	%Önceki hata kovaryansı
	%Dt ynai örnekleme süresi
	%params yani sensör gürültü ve sistem parametreleri
	%Çıktılar
	%Güncel quaternion tahmini
	%güncel gyro bias tahmini
	%güncel hata kovaryansı

	
	
	
	%gyro functionundan gelen ölçümü kalibre ediyoruz
	gyro_body = params.gyro.misalignment' * current_gyro;
	omega_estimated = gyro_body - bias_prev;
	%burada kinematik hesaplanıyor ve çıktı olarak tahmin quaternion 
	[q_pred, ~] = kinematics(q_prev, omega_estimated, dt);


	%açısal hızın skew symmetric matrisi hesaplanıyor
	omega_skew = [0, -omega_estimated(3), omega_estimated(2); ...
		omega_estimated(3), 0, -omega_estimated(1); ...
		-omega_estimated(2), omega_estimated(1), 0];

	%HATA DURUM MODELİ OLUŞTURULUYOR
	F = [-omega_skew, -eye(3); zeros(3), zeros(3)];
	Phi = eye(6) + F * dt;

	%TOPLAM ÖLÇÜM GÜRÜLTÜSÜ OLUŞTURULUYOR
	Q_attitude = params.gyro.covariance * dt;
	Q_bias = params.gyro.bias_random_walk^2 * dt * eye(3);
	Q = blkdiag(Q_attitude, Q_bias);

	%HATA KOVARTANSI TAHMİN EDİLİYOR
	P_pred = Phi * P_prev * Phi' + Q;

	% YILDIZ İZLEYİCİNİN DCMDEN QUATERNİONA DÖNÜŞÜM.
	q_measurement = dcm_to_quaternion(current_star_tracker_dcm);

	%QUATERNYON İŞARET BELİRSİZLİ GİDERME
	if dot(q_measurement, q_pred) < 0
		q_measurement = -q_measurement;
	end

	%TAHMİN İLE ÖLÇÜM ARASINDAKİ QUATERNİYON FARKI HESAPLANIYORT
	q_error = quaternion_multiply(q_measurement, quaternion_conjugate(q_pred));
	if q_error(1) < 0
		q_error = -q_error;
	end
	innovation = 2 * q_error(2:4);

	%ÖLÇÜM MODELİ
	H = [eye(3), zeros(3)];
	R = params.startracker.covariance;

	%İNNOVATİON KOVARYANSI
	S = H * P_pred * H' + R;


	%KALMAN KAZANCI

	K = P_pred * H' / S;
	%HATA VEKTÖRÜNÜN HESAPLANMASI

	correction = K * innovation;

	%KÜÇÜK AÇI DÜZELTMESİ
	delta_theta = correction(1:3);
	%KÜÇÜK AÇI DÜZELTMESİNİN QUATERNİYONA DÖNÜŞMESİ
	delta_q = [1; 0.5 * delta_theta];
	delta_q = delta_q / norm(delta_q);

	%TUTUM DÜZELTMESİ
	q_hat = quaternion_multiply(delta_q, q_pred);
	q_hat = q_hat / norm(q_hat);
%BİAS GÜNCELLEMESİ
	bias_hat = bias_prev + correction(4:6);
	%KOVARYANS GÜNCELLENMESİ
	P_hat = (eye(6) - K * H) * P_pred;
	P_hat = 0.5 * (P_hat + P_hat');
end

%QUATERNİON CARPIMI
function q_product = quaternion_multiply(q_left, q_right)
	product_scalar = q_left(1) * q_right(1) - ...
		dot(q_left(2:4), q_right(2:4));
	product_vector = q_left(1) * q_right(2:4) + ...
		q_right(1) * q_left(2:4) + ...
		cross(q_left(2:4), q_right(2:4));
	q_product = [product_scalar; product_vector(:)];
end

function q_conjugate = quaternion_conjugate(q)
	q_conjugate = [q(1); -q(2:4)];
end
